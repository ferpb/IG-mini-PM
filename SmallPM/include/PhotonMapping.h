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

#ifndef _PHOTONMAPPING_H_
#define _PHOTONMAPPING_H_

#include <time.h>

#include "globals.h"
#include "LightSource.h"
#include "Vector3.h"
#include "KDTree.h"
#include "World.h"
#include "ParticipatingMedia.h"

class World;
class Intersection;
class Ray;

enum KERNEL
{
    NORMAL,
    GAUSSIAN,
    CONE
};
/** The PhotonMapping class is a light integrator that implements
	classic Jensen's photon mapping algorithm. Both the pre-process
	pass and the gathering pass are assumed to be iterative (as 
	oposed to the traditional recursive formulation of Whitted's 
	ray tracing. */

// Structure defining a photon (a directionally-resolved packet of
// energy), that will be used later for radiance estimation.
struct Photon
{
    Vector3 position;
    Vector3 direction;
    Vector3 flux;
    // Se guarda la luz a la que pertenece el photon
    int num_light;

    Photon() : position(0), direction(0), flux(0) {}
    Photon(const Vector3 &p, const Vector3 &d, const Vector3 &f, int _num_light) : position(p), direction(d), flux(f), num_light(_num_light) {}
};

class PhotonMapping
{
    World *world;
    KERNEL kernel;

    unsigned int m_max_nb_shots, m_nb_global_photons, m_nb_caustic_photons, m_nb_volumen_photons;
    unsigned int m_nb_current_shots;

    unsigned int m_nb_photons;
    bool m_raytraced_direct, foggy;
    std::vector<Real> lights_probs;
    std::vector<int> shots_per_light;

    KDTree<Photon, 3> m_global_map, m_caustics_map, m_volumen_photones;

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
    bool trace_ray(const Ray &r, const Vector3 &p,
                   std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, std::list<Photon> &volumen_photones, bool direct, int num_light);

  public:
    PhotonMapping(World *_world, unsigned int nb_global_photons, unsigned int nb_caustic_photons, unsigned int photons_volume,
                  unsigned int max_nb_shots, unsigned int nb_photons, bool raytraced_direct = false, KERNEL _kernel = NORMAL, bool _foggy = false) : world(_world), m_max_nb_shots(max_nb_shots), m_nb_current_shots(0),
                                                                                                                                                     m_nb_global_photons(nb_global_photons), m_nb_caustic_photons(nb_caustic_photons),
                                                                                                                                                     m_nb_photons(nb_photons), m_raytraced_direct(raytraced_direct)
    {

        foggy = _foggy;
        m_nb_volumen_photons = photons_volume;
        kernel = _kernel;
        Real acum = 0.0f;
        int num_lights = world->nb_lights();

        for (int i = 0; i < num_lights; i++)
        {
            acum += world->light(i).get_intensities().avg();
            lights_probs.push_back(world->light(i).get_intensities().avg());
        }

        lights_probs[0] = lights_probs[0] / acum;

        for (int i = 1; i < lights_probs.size(); i++)
            lights_probs[i] = lights_probs[i - 1] + lights_probs[i] / acum;

        shots_per_light = std::vector<int>(num_lights);
    }

    // Preprocess the photon map. This needs to be run before rendering,
    // or no photons will be stored to compute radiance in the rendering
    // pass.
    void preprocess();

    // Computes shading at the intersection 'it0' and returns the estimated
    // radiance.
    Vector3 shade(Ray ray, Intersection &it0) const;
    int Russian_roulette() const;
    Vector3 Estimacion(Vector3 interesec_point, const KDTree<Photon, 3> &photon_map, int nb_photons, std::vector<int> photons_per_light) const;
    Vector3 Direct_light_NEE(const Intersection &it) const;
    Vector3 Direc_light_RR(const Intersection &it) const;
};

vector<Real> vector3_to_vector(const Vector3 &v);

#endif
