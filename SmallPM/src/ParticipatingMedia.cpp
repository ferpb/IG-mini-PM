#include "ParticipatingMedia.h"
#include <cmath>

bool fog_trace(Ray &photon_ray, Vector3 &energy, Intersection &it)
{
    /*
    Real distance = random_real() * photon_ray.get_parameter();

    if (1 - exp(-distance * sigma_t) > random_real())
    {

        Ray out_ray = Ray(photon_ray.get_origin(), photon_ray.get_direction(), photon_ray.get_level());
        out_ray.cond_set_parameter(distance);
        it = Intersection(out_ray, new Sphere({}, 1, new Fog(uniform_sphere_sample(), sigma_s / sigma_t)), Vector3(), Vector2());

        return true;
    }*/

    Real distance = -1 * (log(random_real(1e-10, 1.0 - 1e-10)) / sigma_t);
    if (1 - exp(-distance * sigma_t) > random_real())
    {
        Ray out_ray = Ray(photon_ray.get_origin(), photon_ray.get_direction(), photon_ray.get_level());
        out_ray.cond_set_parameter(distance);
        it = Intersection(out_ray, new Sphere({}, 1, new Fog(uniform_sphere_sample(), sigma_s / sigma_t)), Vector3(), Vector2());

        return true;
    }

    else
    {
        return false;
    }
}

Vector3 fog_ray_marching(const Vector3 &direction, const Vector3 &origin, const Real ray_length, int nb_volumen_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen)
{

    Vector3 Ld, Li;
    Vector3 x = origin;

    Real length_partition = 0.07f;
    Real acum_length = length_partition;
    Real radius;
    std::vector<const KDTree<Photon, 3>::Node *> nodes;

    if (!m_volumen_photones.is_empty())
    {
        while (acum_length < ray_length)
        {

            x = x + (direction * length_partition);
            m_volumen_photones.find(vector3_to_vector(x), 50, nodes, radius);
            for (const auto &node : nodes)
            {
                Li += node->data().flux / (float)shots_per_light_volumen[node->data().num_light];
            }
            Li *= (1.0f / (4.0f * M_PI)) * (3.0f / (4.0f * M_PI * radius * radius * radius));
            Li *= (sigma_s / sigma_t);

            Ld = Ld * exp(-length_partition * sigma_t) + sigma_s * length_partition * Li;
            acum_length += length_partition;
            Li = Vector3(0, 0, 0);
        }
    }
    return Ld;
}

Vector3 fog_ray_marching_probs(const Vector3 &direction, const Vector3 &origin, const Real ray_length, int nb_volumen_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen)
{

    Vector3 Ld, Li;
    Vector3 x = origin;
    Vector3 x_prim;

    Real length_partition = ray_length / num_partitions;
    Real acum_length = length_partition;
    Real radius;
    Real length_t;
    std::vector<const KDTree<Photon, 3>::Node *> nodes;

    if (!m_volumen_photones.is_empty())
    {
        while (acum_length < ray_length)
        {

            x_prim = x + (direction * (length_partition * random_real(0, 1.0f - 1e-100)));
            m_volumen_photones.find(vector3_to_vector(x_prim), nb_volumen_photons, nodes, radius);
            for (const auto &node : nodes)
            {
                Li += node->data().flux / (float)shots_per_light_volumen[node->data().num_light];
            }
            Li *= (1.0f / (4.0f * M_PI)) * (3.0f / (4.0f * M_PI * radius * radius * radius));
            Li *= (sigma_s / sigma_t);
            length_t = (x_prim - x).length();
            //Ld = Ld * exp(-length_t * sigma_t) + sigma_s * length_t * Li;

            Ld = Ld * exp(-length_t * sigma_t) + sigma_s * Li * length_t;
            acum_length += length_partition;
            Li = Vector3(0, 0, 0);
            x = x_prim;
        }
    }
    return Ld;
}

Vector3 fog_ray_marching_probs2(const Vector3 &direction, const Vector3 &origin, const Real ray_length, int nb_volumen_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen)
{

    Vector3 Ld, Li;
    Vector3 x = origin;
    Vector3 x_prim;

    Real length_t = -1 * (log(random_real(1e-10, 1.0 - 1e-10)) / sigma_t);
    Real acum_length = length_t;
    Real radius;

    std::vector<const KDTree<Photon, 3>::Node *> nodes;

    if (!m_volumen_photones.is_empty())
    {
        while (acum_length < ray_length)
        {

            x_prim = x + (direction * length_t);
            m_volumen_photones.find(vector3_to_vector(x_prim), nb_volumen_photons, nodes, radius);
            for (const auto &node : nodes)
            {
                Li += node->data().flux / (float)shots_per_light_volumen[node->data().num_light];
            }
            Li *= (1.0f / (4.0f * M_PI)) * (3.0f / (4.0f * M_PI * radius * radius * radius));
            Li *= (sigma_s / sigma_t);

            //Ld = Ld * exp(-length_t * sigma_t) + sigma_s * length_t * Li;

            Ld = Ld * exp(-length_t * sigma_t) + sigma_s * Li * length_t;
            length_t = -1 * (log(random_real(1e-10, 1.0 - 1e-10)) / sigma_t);
            acum_length += length_t;
            Li = Vector3(0, 0, 0);
            x = x_prim;
        }
    }
    return Ld;
}
Vector3 attenuation(Real distance, Vector3 Lum)
{
    return exp(-distance * sigma_t) * Lum;
}
