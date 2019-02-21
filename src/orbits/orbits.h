#ifndef ORBITS_H
#define ORBITS_H

#include "../type_class.h"
#include "../own_math.h"
#include "../macros.h"
#include <math.h>

namespace SpaceH::obt {

    template<typename Scalar>
    inline Scalar myacos(Scalar x) {
        return acos(SpaceH::min(SpaceH::max(-1.0, x), 1.0));
    }

    template<typename Vector, typename Scalar>
    Vector calcu_eccentricity(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                              const Vector &vel2) {
        Vector dr = pos1 - pos2;
        Vector dv = vel1 - vel2;
        Scalar u = Const::G * (m1 + m2);
        return (dr * (dv.norm2() - u * dr.reNorm()) - dv * dot(dr, dv)) / u;
    }


    template<typename Vector, typename Scalar>
    Scalar calcu_semi_major_axis(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                                 const Vector &vel2) {
        Vector dr = pos1 - pos2;
        Vector dv = vel1 - vel2;
        Scalar u = Const::G * (m1 + m2);
        Scalar r = dr.norm();
        return -u * r / (r * dv.norm2() - 2 * u);
    }

    template<typename Vector, typename Scalar>
    void euler_rotate(Vector &v, const Scalar phi, const Scalar theta, const Scalar psi) {
        Scalar sin_phi = sin(phi);
        Scalar cos_phi = cos(phi);
        Scalar sin_psi = sin(psi);
        Scalar cos_psi = cos(psi);
        Scalar sin_theta = sin(theta);
        Scalar cos_theta = cos(theta);

        Scalar x = v.x * (cos_phi * cos_psi - sin_phi * cos_theta * sin_psi)
                   - v.y * (cos_phi * sin_psi + sin_phi * cos_theta * cos_psi)
                   + v.z * (sin_phi * sin_theta);

        Scalar y = v.x * (sin_phi * cos_psi + cos_phi * cos_theta * sin_psi)
                   - v.y * (sin_phi * sin_psi - cos_phi * cos_theta * cos_psi)
                   - v.z * (cos_phi * sin_theta);

        Scalar z = v.x * sin_theta * sin_psi + v.y * sin_theta * cos_psi + v.z * cos_theta;

        v.x = x, v.y = y, v.z = z;
    };

    template<typename Scalar>
    Scalar get_random_mean_anomaly(Scalar e, Scalar Mmin, Scalar Mmax) {
        if (e >= 0)
            return Random<Scalar>::uniform() * (Mmax - Mmin) + Mmin;
        else {
            SPACEHUB_ABORT("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }

    template<typename Scalar>
    Scalar calcu_true_anomaly(Scalar E, Scalar e) {
        if (0 <= e && e < 1)
            return 2 * atan2(sqrt(1 + e) * sin(E * 0.5), sqrt(1 - e) * cos(0.5 * E));
        else if (e > 1)
            return 2 * atan2(sqrt(1 + e) * sinh(E * 0.5), sqrt(e - 1) * cosh(0.5 * E));
        else if (iseq(e, 1.0))
            return 2 * atan(0.5 * E);
        else {
            SPACEHUB_ABORT("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }


    template<typename Scalar>
    Scalar calcu_eccentric_anomaly(Scalar M, Scalar e) {
        if (0 <= e && e < 1)//newton iteration may encounter stationary point
            return SpaceH::root_dichotom(
                    [&](Scalar x) -> Scalar { return x - e * sin(x) - M; });//find this function in ownMath.h
        else if (e > 1)
            return SpaceH::root_dichotom([&](Scalar x) -> Scalar { return e * sinh(x) - x - M; });
        else if (fabs(e - 1) < SpaceH::epsilon<Scalar>::value)
            return SpaceH::root_dichotom([&](Scalar x) -> Scalar { return x + x * x * x / 3 - M; });
        else {
            SPACEHUB_ABORT("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }

    enum class OrbitType {
        ellipse, parabola, hyperbola, none
    };

    template<typename T>
    OrbitType classify_orbit(T eccentricity) {
        if (0 <= eccentricity && eccentricity < 1) {
            return OrbitType::ellipse;
        } else if (iseq(eccentricity, 1.0)) {
            return OrbitType::parabola;
        } else if (eccentricity > 1) {
            return OrbitType::hyperbola;
        } else {
            return OrbitType::none;
        }
    }

    struct RandomIndicator {
    } thermal;

    template<typename Real>
    struct OrbitParam {
        using Scalar = Real;
        Scalar m1;
        Scalar m2;
        Scalar e;//eccentricity
        Scalar p;//semi-latus rectum
        Scalar i;//inclination
        Scalar Omega;//longitude of the ascending node
        Scalar omega;//argument of periapsis
        Scalar nu;//true anomaly
        OrbitType orbit_type;

        OrbitParam() = delete;

        template<typename Variant1, typename Variant2, typename Variant3, typename Variant4>
        OrbitParam(Scalar _m1_, Scalar _m2_, Scalar _p_, Scalar _e_, Variant1 _inclination_, Variant2 _phi_,
                   Variant3 _psi_, Variant4 _trueAnomaly_) {
            if (_p_ < 0) SPACEHUB_ABORT("Semi-latus rectum cannot be negative");

            orbit_type = classify_orbit(_e_);

            if (orbit_type == OrbitType::none) SPACEHUB_ABORT("Eccentrcity cannot be negative or NaN!");

            m1 = _m1_;
            m2 = _m2_;
            p = _p_;
            e = _e_;

            if constexpr (std::is_same_v<Variant1, Scalar>) {
                i = _inclination_;
            } else if (std::is_same_v<Variant1, RandomIndicator>) {
                shuffle_i();
            } else {
                SPACEHUB_ABORT("Unexcepted type of inclination!");
            }

            if constexpr (std::is_same_v<Variant2, Scalar>) {
                Omega = _phi_;
            } else if (std::is_same_v<Variant2, RandomIndicator>) {
                shuffle_Omega();
            } else {
                SPACEHUB_ABORT("Unexcepted type of longitude of the ascending node!");
            }

            if constexpr (std::is_same_v<Variant3, Scalar>) {
                omega = _psi_;
            } else if (std::is_same_v<Variant3, RandomIndicator>) {
                shuffle_omega();
            } else {
                SPACEHUB_ABORT("Unexcepted type of argument of periapsis!");
            }

            if constexpr (std::is_same_v<Variant4, Scalar>) {
                nu = _trueAnomaly_;
            } else if (std::is_same_v<Variant4, RandomIndicator>) {
                shuffle_nu();
            } else {
                SPACEHUB_ABORT("Unexcepted type of true_anomaly!");
            }
        }

        inline void shuffle_i() {
            i = acos(Random<Scalar>::uniform() * 2 - 1);
        }

        inline void shuffle_Omega() {
            Omega = Random<Scalar>::uniform() * 2 * Const::PI - Const::PI;
        }

        inline void shuffle_omega() {
            omega = Random<Scalar>::uniform() * 2 * Const::PI - Const::PI;
        }

        inline void shuffle_nu() {
            if (orbit_type == OrbitType::ellipse) {
                Scalar M = obt::get_random_mean_anomaly(e, -Const::PI, Const::PI);
                Scalar E = obt::calcu_eccentric_anomaly(M, e);
                nu = obt::calcu_true_anomaly(E, e);
            } else {
                SPACEHUB_ABORT("Only elliptical orbit provides random anomaly method at this moment!");
            }
        }

        friend std::ostream &operator<<(std::ostream &os, OrbitParam const &obt) {
            SpaceH::display(os, obt.m1, obt.m2, obt.e, obt.p, obt.i, obt.Omega, obt.omega, obt.nu);
            return os;
        }
    };

    template<typename Particle, typename Scalar>
    void place_at_orbit(Particle &particle, OrbitParam<Scalar> const &param) {
        Scalar u = (param.m1 + param.m2) * Const::G;

        Scalar sin_nu = sin(param.nu);
        Scalar cos_nu = cos(param.nu);

        Scalar r = param.p / (1 + param.e * cos_nu);
        Scalar v = sqrt(u / param.p);

        particle.pos = r * Vec3<Scalar>(cos_nu, sin_nu, 0);
        particle.vel = v * Vec3<Scalar>(-sin_nu, param.e + cos_nu, 0);

        obt::euler_rotate(particle.pos, param.Omega, param.i, param.omega + Const::PI);
        obt::euler_rotate(particle.vel, param.Omega, param.i, param.omega + Const::PI);
    }

    template<typename Particle>
    struct Kepler {
        /* Typedef */
        using Scalar = typename Particle::Scalar;
        using Vector = typename Particle::Vector;
        using ParticleType = Particle;

        /* Typedef */

        Kepler() = delete;

        Kepler(OrbitParam<Scalar> const &param) {
            if (param.m1 >= param.m2) {
                obj1_.mass = param.m1;
                obj2_.mass = param.m2;
            } else {
                obj1_.mass = param.m2;
                obj2_.mass = param.m1;
            }

            obj1_.pos = Vector(0.0, 0.0, 0.0);
            obj1_.vel = Vector(0.0, 0.0, 0.0);

            Scalar u = (param.m1 + param.m2) * Const::G;

            Scalar sin_nu = sin(param.nu);
            Scalar cos_nu = cos(param.nu);

            Scalar r = param.p / (1 + param.e * cos_nu);
            Scalar v = sqrt(u / param.p);

            obj2_.pos = r * Vector(cos_nu, sin_nu, 0);
            obj2_.vel = v * Vector(-sin_nu, param.e + cos_nu, 0);

            obt::euler_rotate(obj2_.pos, param.Omega, param.i, param.omega + Const::PI);
            obt::euler_rotate(obj2_.vel, param.Omega, param.i, param.omega + Const::PI);
        }

        template<typename Variant1, typename Variant2, typename Variant3, typename Variant4>
        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e, Variant1 i , Variant2 phi, Variant3 psi, Variant4 trueAnomaly)
            : Kepler(OrbitParam<Scalar>(m1, m2, p, e, i, phi, psi, trueAnomaly)) {}

        Particle &first() {
            return obj1_;
        }

        Particle &second() {
            return obj2_;
        }

        void move_CoM_to(Particle const &partc) {
            move_to_CM_coord();
            obj1_.pos += partc.pos;
            obj1_.vel += partc.vel;
            obj2_.pos += partc.pos;
            obj2_.vel += partc.vel;
        }

        void move_CoM_to(Vector const &CMpos, Vector const &CMvel) {
            move_to_CM_coord();
            obj1_.pos += CMpos;
            obj1_.vel += CMvel;
            obj2_.pos += CMpos;
            obj2_.vel += CMvel;
        }

        void move_to_CM_coord() {
            Scalar tot_mass = obj1_.mass + obj2_.mass;
            Vector CMpos((obj1_.pos * obj1_.mass + obj2_.pos * obj2_.mass) / tot_mass);
            Vector CMvel((obj1_.vel * obj1_.mass + obj2_.vel * obj2_.mass) / tot_mass);
            obj1_.pos -= CMpos;
            obj1_.vel -= CMvel;
            obj2_.pos -= CMpos;
            obj2_.vel -= CMvel;
        }

    private:
        Particle obj1_;
        Particle obj2_;
    };


    template<typename Container>
    void create_particle_id(Container &particles) {
        size_t num = particles.size();
        size_t i = 0;
        for (auto &p : particles) {
            p.idn = i++;
        }
    }

    template<typename Vector, typename Scalar>
    void to_oribt_parameters(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                             const Vector &vel2, OrbitParam<Scalar> &params) {
        Vector dr = pos1 - pos2;
        Vector dv = vel1 - vel2;
        Vector L = cross(dr, dv);
        Vector N = cross(Vector(0, 0, 1.0), L);
        Scalar u = Const::G * (m1 + m2);
        Scalar r = dr.norm();
        Scalar n = N.norm();
        Scalar l = L.norm();
        Scalar rv = dot(dr, dv);
        Vector E = (dr * (dv.norm2() - u * dr.reNorm()) - dv * rv) / u;
        params.e = E.norm();

        if (!iseq(params.e, 1.0)) {
            Scalar a = -u * r / (r * dv.norm2() - 2.0 * u);
            params.p = a * (1.0 - params.e * params.e);
        } else {
            params.p = l * l / u;
        }

        params.i = acos(L.z / l);

        if (params.e != 0) {
            params.nu = SpaceH::sign(rv) * myacos(dot(E / params.e, dr / r));

            if (n != 0) {
                params.Omega = SpaceH::sign(N.y) * myacos(N.x / n);
                params.omega = SpaceH::sign(E.z) * myacos(dot(E / params.e, N / n));
            } else {
                params.omega = -SpaceH::sign(E.y) * myacos(-E.x / params.e);
                params.Omega = params.omega;
            }
        } else {
            if (n != 0) {
                params.Omega = SpaceH::sign(N.y) * myacos(N.x / n);
                params.omega = 0;
                Vector peri = cross(L, N);
                params.nu = -SpaceH::sign(dot(N, dr)) * myacos(dot(peri / peri.norm(), dr / r));
            } else {
                params.Omega = params.omega = 0;
                params.nu = SpaceH::sign(dr.y) * acos(dot(Vector(1.0, 0, 0), dr / r));
            }
        }
    }
}
#endif
