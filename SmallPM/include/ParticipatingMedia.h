#pragma once

#include "KDTree.h"
#include "Vector3.h"
#include <cmath>
#include "Ray.h"
#include "Intersection.h"
#include "Sphere.h"
#include "Random.h"
#include "BSDF.h"
#include "PhotonMapping.h"
#include <iostream>
#include <iterator>
#include "World.h"

#define sigma_t 8.0f
#define sigma_s 0.5f
#define num_partitions 15.0f

struct Photon;

class Fog : public BSDF
{

    Vector3 next_ray;
    float albedo;

  public:
    Fog(const Vector3 _next_ray, float _albedo) : BSDF(nullptr), next_ray(_next_ray), albedo(_albedo){};

    void get_outgoing_sample_ray(const Intersection &it, Ray &r, Real &pdf) const
    {

        r = Ray(it.get_position(), this->next_ray, 0);
        pdf = 1;
    };
    Vector3 get_albedo(const Intersection &it) const
    {

        return Vector3(albedo, albedo, albedo);
    };
    Real get_specular(const Intersection &it) const
    {

        return 0.0f;
    };
    bool is_delta() const
    {

        return false;
    };
};

bool fog_trace(Ray &photon_ray, Vector3 &energy, Intersection &it);
Vector3 fog_ray_marching(const Vector3 &direction, const Vector3 &origin, const Real distance, int m_nb_caustic_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen);
Vector3 fog_ray_marching_probs(const Vector3 &direction, const Vector3 &origin, const Real distance, int m_nb_caustic_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen);
Vector3 fog_ray_marching_probs2(const Vector3 &direction, const Vector3 &origin, const Real ray_length, int nb_volumen_photons, const KDTree<Photon, 3> &m_volumen_photones, std::vector<int> shots_per_light_volumen);
Vector3 attenuation(Real distance, Vector3 Lum);
