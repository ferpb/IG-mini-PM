/*********************************************************************************
 * Random number generation functions 
 *
 * File: random.cpp
 * Author: Fernando Peña (NIA: 756012)
 * Author: Jose Daniel Subias Sarrato (NIA: 759533)
 * Date: 1/12/2020
 * Coms: Informática Gráfica, 2020-2021
 **********************************************************************************/

#include "Random.h"

#include <random>
#include <cmath>

static thread_local std::random_device rd; // obtain a random number from hardware
static thread_local std::mt19937 gen(rd()); // seed the generator
static thread_local std::uniform_real_distribution<Real> distr(0.0f, 1.0f); // define the range

Real random_real()
{
    return distr(gen);
}

/* Scaling a number x in the range [a, b] to [min, max]:
 *
 *     f(x) = ((max-min)(x-a) / (b-a)) + min
 *
 */

Real random_real(Real n)
{
    return n * distr(gen);
}

Real random_real(Real min, Real max)
{
    return (max - min) * distr(gen) + min;
}

int random_int()
{
    return distr(gen);
}

int random_int(int n)
{
    return n * distr(gen);
}

int random_int(int min, int max)
{
    return (max - min) * distr(gen) + min;
}

Vector3 uniform_sphere_sample()
{
    // Real inclination = M_PI * random_real(); // theta
    // Real inclination = acosf(2 * random_real() - 1); // theta
    // Real inclination = acosf(sqrtf(random_real())); // theta
    // Real a = -1.0f;
    // Real b = 1.0f;

    Real inclination = acosf(random_real(-1, 1)); // theta
    //Real inclination = acosf(random_real()); // theta
    Real azimuth = 2 * M_PI * random_real(); // phi

    // return Vector3(sinf(inclination) * sinf(azimuth), sinf(inclination) * cosf(azimuth), cosf(inclination));
    return Vector3(sinf(inclination) * cosf(azimuth), sinf(inclination) * sinf(azimuth), cosf(inclination));
}

/*
#define CDF
// #define ARCHIMEDES
// #define REJECTION
// #define REJECTION_JENSEN

#if defined(CDF)
*/

/*
#elif defined(ARCHIMEDES)
Vector3 uniform_sphere_sample()
{
    // Utilizando el teorema de arquímedes
    Real theta = 2 * M_PI * random_real();
    Real z = random_real();
    Real x = sqrtf(1 - z * z) * cosf(theta);
    Real y = sqrtf(1 - z * z) * sinf(theta);
    return Vector3(x, y, z);
}

#elif defined(REJECTION)
Vector3 uniform_sphere_sample()
{
    Real x, y, z;
    do
    {
        x = random_real(-1, 1);
        y = random_real(-1, 1);
        z = random_real(-1, 1);
        // } while (pow(x, 2) + pow(y, 2) + pow(z, 2) > 1);
    } while (x * x + y * y + z * z > 1);
    return Vector3(x, y, z);
}

#elif defined(REJECTION_JENSEN)
Vector3 uniform_sphere_sample()
{
    float x, y, z;
    do
    {
        // Esto es lo mismo que hace random_real sobrecargada
        x = 2 * random_real() - 1;
        y = 2 * random_real() - 1;
        z = 2 * random_real() - 1;
    } while (x * x + y * y + z * z > 1);

    return Vector3(x, y, z);
}

#endif
*/
/*
Vector3 uniform_sphere_sample()
{
    Real x, y, z;
    do
    {
        x = random_real(-1, 1);
        y = random_real(-1, 1);
        z = random_real(-1, 1);
        // } while (pow(x, 2) + pow(y, 2) + pow(z, 2) > 1);
    } while (x * x + y * y + z * z > 1);
    return Vector3(x, y, z);
}
*/
