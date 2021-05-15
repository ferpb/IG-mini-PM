/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not 
intend to be fast or general, but just to provide an educational tool for undergraduate
students. 

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"
#include "BVH.h"

#include "Random.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>

vector<Real> vector3_to_vector(const Vector3 &v)
{
    return vector<Real>{v.getComponent(0), v.getComponent(1), v.getComponent(2)};
}

inline std::ostream &operator<<(std::ostream &os, const Vector3 &v)
{
    os << "[ " << v.getComponent(0) << ", " << v.getComponent(1) << ", " << v.getComponent(2) << " ]";
    return os;
}

#define EPSILON 0.01f

bool similar(Real a, Real b)
{
    return fabs(a - b) <= EPSILON;
}

bool similar(Vector3 a, Vector3 b)
{
    return similar(a.getComponent(0), b.getComponent(0)) && similar(a.getComponent(1), b.getComponent(1)) && similar(a.getComponent(2), b.getComponent(2));
}

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct
// photons (from light to surface) are being stored or not.
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray &r, const Vector3 &p,
                              std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, std::list<Photon> &volumen_photones, bool direct, int num_light)
{

    //Check if max number of shots done...
    if (++m_nb_current_shots > m_max_nb_shots)
    {
        return false;
    }

    // Compute irradiance photon's energy
    Vector3 energy(p);

    Ray photon_ray(r);
    photon_ray.shift();

    bool is_caustic_particle = false;
    bool is_foggy = false;
    int initial = 0;

    shots_per_light[num_light] = shots_per_light[num_light] + 1;

    //Iterate the path
    while (1)
    {
        // Throw ray and update current_it
        Intersection it;
        world->first_intersection(photon_ray, it);

        if (!it.did_hit())
            break;

        if (foggy)
            is_foggy = fog_trace(photon_ray, energy, it);

        //Check if has hit a delta material...
        if (it.intersected()->material()->is_delta())
        {
            // If delta material, then is caustic...
            // Don't store the photon!
            is_caustic_particle = true;
        }
        else if (photon_ray.get_level() > 0 || direct || is_foggy)
        {
            //If non-delta material, store the photon!

            if (is_foggy)
            {

                if (volumen_photones.size() < m_nb_volumen_photons)
                    volumen_photones.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy, num_light));
            }
            else if (is_caustic_particle)
            {

                //If caustic particle, store in caustics
                if (caustic_photons.size() < m_nb_caustic_photons)
                    caustic_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy, num_light));
            }
            else
            {

                //If non-caustic particle, store in global
                if (global_photons.size() < m_nb_global_photons)
                    global_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy, num_light));
            }
            is_caustic_particle = false;
        }

        Real pdf;

        Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
        Real avg_surf_albedo = surf_albedo.avg();

        Real epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
        while (epsilon2 < 0.)
            epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

        if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20)
            break;

        // Random walk's next step
        // Get sampled direction plus pdf, and update attenuation
        it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);

        // Shade...
        energy = energy * surf_albedo;
        if (!it.intersected()->material()->is_delta())
        {
            if (!is_foggy)
                energy *= dot_abs(it.get_normal(), photon_ray.get_direction()) / 3.14159;
        }

        energy = energy / (pdf * avg_surf_albedo);
    }

    if (caustic_photons.size() == m_nb_caustic_photons &&
        global_photons.size() == m_nb_global_photons && volumen_photones.size() == m_nb_volumen_photons)
    {
        m_max_nb_shots = m_nb_current_shots - 1;
        return false;
    }

    return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering.
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{

    // Stores true if there are more photons to be shot

    std::list<Photon> global_photons = {};
    std::list<Photon> caustic_photons = {};
    std::list<Photon> volumen_photones = {};

    int idx;
    const LightSource *ls;
    Vector3 dir;

    do
    {

        // Chose a random light source
        idx = Russian_roulette();
        //std::cout << "rr: " << idx << std::endl;
        ls = &world->light(idx);

        // Sample a random point on the unit sphere
        dir = uniform_sphere_sample();

        // Shot photons while there are photons to be shot
    } while (trace_ray(Ray(ls->get_position(), dir.normalize()), (ls->get_intensities() * (4.0f * M_PI)), global_photons, caustic_photons, volumen_photones, !m_raytraced_direct, idx));

    // Creo que get intesities si que hay que dividirlo, porque esta es la energía total de la fuente de luz,
    // luego cada fotón tendrá una posción de esa energía.
    // Pero no se si este cálculo bien

    /*
    plot_photons(global_photons, "maps/global");
    plot_photons(caustic_photons, "maps/caustic");
    plot_photons(volumen_photones, "maps/volume");

    std::list<Photon> global_caustic_photons = global_photons;
    global_caustic_photons.insert(global_caustic_photons.end(), caustic_photons.begin(), caustic_photons.end());
    plot_photons(global_caustic_photons, "maps/global_caustic");

    std::list<Photon> all_photons = global_photons;
    all_photons.insert(all_photons.end(), caustic_photons.begin(), caustic_photons.end());
    all_photons.insert(all_photons.end(), volumen_photones.begin(), volumen_photones.end());
    plot_photons(all_photons, "maps/all");
    */

    // Store the photons in the kd-trees
    m_global_map.clear();
    m_caustics_map.clear();
    m_volumen_photones.clear();

    if (!global_photons.empty())
    {
        for (const Photon &photon : global_photons)
        {
            m_global_map.store(vector3_to_vector(photon.position), photon);
        }
        m_global_map.balance();
        std::cout << "global: " << m_global_map.nb_elements() << std::endl;
    }

    if (!caustic_photons.empty())
    {

        for (const Photon &photon : caustic_photons)
        {
            m_caustics_map.store(vector3_to_vector(photon.position), photon);
        }
        m_caustics_map.balance();
        std::cout << "causticas:  " << m_caustics_map.nb_elements() << std::endl;
    }

    if (!volumen_photones.empty())
    {

        for (const Photon &photon : volumen_photones)
        {
            m_volumen_photones.store(vector3_to_vector(photon.position), photon);
        }
        m_volumen_photones.balance();
        std::cout << "photones en el volumen:  " << m_volumen_photones.nb_elements() << std::endl;
    }

    for (int i = 0; i < shots_per_light.size(); i++)
    {

        std::cout << "nb light " << i << ": " << shots_per_light[i] << std::endl;
    }

    std::cout << "Emission finished" << std::endl;
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Ray ray, Intersection &it0) const
{

    Intersection it(it0);
    Vector3 L_global(0), L_causticas(0), L_direct(0), L_fog(0);
    Real total_distance = ray.get_parameter();

    Vector3 W(1);

    int bounces = 0;
    // If the material is delta, follow the path until we hit a diffuse material
    while (it.intersected()->material()->is_delta() && bounces < 20)
    {

        // pdf is not used in delta materials
        Real pdf;
        it.intersected()->material()->get_outgoing_sample_ray(it, ray, pdf);
        total_distance += ray.get_parameter();
        if (!it.did_hit())
        {
            return Vector3(0);
        }

        ray.shift();
        world->first_intersection(ray, it);

        bounces++;
    }

    if (!it.did_hit() || it.intersected()->material()->is_delta())
    {
        return Vector3(0);
    }

    if (foggy)
        L_fog += fog_ray_marching_probs2(ray.get_direction(), ray.get_origin(), ray.get_parameter(), m_nb_photons, m_volumen_photones, shots_per_light);

    Vector3 it_position = it.get_position();
    Vector3 it_albedo = it.intersected()->material()->get_albedo(it);
    Vector3 it_normal = it.get_normal();

    if (m_raytraced_direct)
    {
        // Perform next event estimation
        L_direct = Direct_light_RR(it);
    }

    // Radiance estimation using the photon maps
    L_global = Estimation(it_position, m_global_map, m_nb_photons, shots_per_light);
    L_global = L_global * (it_albedo / M_PI);
    L_causticas = Estimation(it_position, m_caustics_map, m_nb_photons, shots_per_light);
    L_causticas = L_causticas * (it_albedo / M_PI);

    if (foggy)
    {

        return L_fog + attenuation(total_distance, L_global + L_causticas + L_direct);
    }
    else
    {

        return L_global + L_causticas + L_direct;
    }
}

Vector3 PhotonMapping::Estimation(Vector3 interesec_point, const KDTree<Photon, 3> &photon_map, int nb_photons, std::vector<int> photons_per_light) const
{
    Vector3 L(0);
    Real radius;
    std::vector<const KDTree<Photon, 3>::Node *> nodes;

    if (!photon_map.is_empty())
    {
        photon_map.find(vector3_to_vector(interesec_point), nb_photons, nodes, radius);

        switch (this->kernel)
        {
        case NORMAL:
        {

            for (const auto &node : nodes)
            {

                L += node->data().flux / (float)photons_per_light[node->data().num_light];
            }

            L *= (1.0 / (M_PI * radius * radius));

            break;
        }
        case CONE:
        {

            Real k = 2.0;
            Real wpc, dp;

            for (const auto &node : nodes)
            {
                dp = (node->data().position - interesec_point).length();
                wpc = 1.0f - dp / (radius * k);
                L += (node->data().flux / (float)photons_per_light[node->data().num_light]) * wpc;
            }

            L *= (1.0 / ((M_PI * radius * radius) * (1.0f - (2.0f / (3.0f * k)))));
            break;
        }
        case GAUSSIAN:
        {

            Real alpha = 0.918;
            Real beta = 1.953;
            Real wpg, dp;
            for (const auto &node : nodes)
            {
                dp = (node->data().position - interesec_point).length();
                wpg = alpha * (1.0 - ((1.0 - exp(-beta * ((dp * dp) / 2.0 * radius * radius))) / (1.0 - exp(-beta))));
                L += (node->data().flux / (float)photons_per_light[node->data().num_light]) * wpg;
            }

            L *= (1.0 / (M_PI * radius * radius));
            break;
        }
        }
    }

    return L;
}

int PhotonMapping::Russian_roulette() const
{
    Real rand_num = random_real();

    for (int i = 0; i < lights_probs.size(); i++)
    {
        if (rand_num <= lights_probs[i])
        {
            return i;
        }
    }

    return -1;
};

Vector3 PhotonMapping::Direct_light_NEE(const Intersection &it) const
{
    Vector3 L_direct(0), L_direct_aux(0);
    Vector3 it_position = it.get_position();
    Vector3 it_albedo = it.intersected()->material()->get_albedo(it);
    Vector3 it_normal = it.get_normal();

    for (int i = 0; i < world->nb_lights(); i++)
    {
        const LightSource *ls = &world->light(i);
        if (ls->is_visible(it_position))
        {

            L_direct_aux += ((it_albedo / M_PI) * ls->get_incoming_light(it_position) * dot_abs(it_normal, (ls->get_position() - it_position).normalize()));
            if (foggy)
            {

                Vector3 shadow_ray = ls->get_position() - it_position;
                Real lenght = shadow_ray.length();
                shadow_ray.normalize();
                Vector3 L_fog_dir_1 = fog_ray_marching_probs2(shadow_ray, ls->get_position(), lenght, m_nb_photons, m_volumen_photones, shots_per_light);
                Vector3 L_fog_dir_2 = attenuation(lenght, L_direct);
                L_direct_aux = L_fog_dir_1 + L_fog_dir_2;
            }
        }
        L_direct += L_direct_aux;
    }
    return L_direct;
}

Vector3 PhotonMapping::Direct_light_RR(const Intersection &it) const
{
    Vector3 L_direct(0);
    Vector3 it_position = it.get_position();
    Vector3 it_albedo = it.intersected()->material()->get_albedo(it);
    Vector3 it_normal = it.get_normal();
    const LightSource *ls = &world->light(Russian_roulette());
    if (ls->is_visible(it_position))
    {

        L_direct += ((it_albedo / M_PI) * ls->get_incoming_light(it_position) * dot_abs(it_normal, (ls->get_position() - it_position).normalize()));
        if (foggy)
        {

            Vector3 shadow_ray = ls->get_position() - it_position;
            Real lenght = shadow_ray.length();
            shadow_ray.normalize();
            Vector3 L_fog_dir_1 = fog_ray_marching_probs2(shadow_ray, ls->get_position(), lenght, m_nb_photons, m_volumen_photones, shots_per_light);
            Vector3 L_fog_dir_2 = attenuation(lenght, L_direct);
            L_direct = L_fog_dir_1 + L_fog_dir_2;
        }
    }

    return L_direct;
}

void PhotonMapping::plot_photons(list<Photon> l, std::string file)
{
    std::ofstream of(file + ".csv");
    for (auto photon : l)
    {

        of << photon.position.getComponent(0) << ", "
           << photon.position.getComponent(1) << ", "
           << photon.position.getComponent(2)
           << std::endl;
    }
}
