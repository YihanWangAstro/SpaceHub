#ifndef ORBITS_H
#define ORBITS_H

#include "../rand-generator.tpp"
#include "../own-math.h"
#include "../macros.h"
#include <math.h>
#include <variant>

namespace SpaceH::Orbit {

    template<typename Scalar>
    inline Scalar myacos(Scalar x) {
        return acos(SpaceH::min(SpaceH::max(-1.0, x), 1.0));
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
    Scalar random_mean_anomaly(Scalar e, Scalar Mmin, Scalar Mmax) {
        if (e >= 0)
            return Random::Uniform<Scalar>::get(Mmin, Mmax);
        else {
            SPACEHUB_ABORT("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }

    template<typename Scalar>
    Scalar calc_true_anomaly(Scalar E, Scalar e) {
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
    Scalar calc_eccentric_anomaly(Scalar M, Scalar e) {
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
    struct OrbitArgs {
    private:
        using Variant = std::variant<Real, RandomIndicator>;
    public:
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

        OrbitArgs() = delete;

        OrbitArgs(Scalar _m1_, Scalar _m2_, Scalar _p_, Scalar _e_, Variant _inclination_, Variant _phi_, Variant _psi_, Variant _true_anomaly_) {
            if (_p_ < 0) SPACEHUB_ABORT("Semi-latus rectum cannot be negative");

            orbit_type = classify_orbit(_e_);

            if (orbit_type == OrbitType::none) SPACEHUB_ABORT("Eccentrcity cannot be negative or NaN!");

            m1 = _m1_;
            m2 = _m2_;
            p = _p_;
            e = _e_;

            auto alt_set = [](Scalar& dst, Scalar& src, auto& method){
                if constexpr (std::holds_alternative<Scalar>(src)) {
                    dst = src;
                } else {
                    method();
                }
            };

            alt_set(i, _inclination_, this->shuffle_i);
            alt_set(Omega, _phi_, this->shuffle_Omega);
            alt_set(omega, _psi_, this->shuffle_omega);
            alt_set(nu, _true_anomaly_, this->shuffle_nu);
        }

        inline void shuffle_i() {
            i = acos(Random::Uniform<Scalar>::get(-1, 1));
        }

        inline void shuffle_Omega() {
            Omega = Random::Uniform<Scalar>::get(-Const::PI, Const::PI);
        }

        inline void shuffle_omega() {
            omega = Random::Uniform<Scalar>::get(-Const::PI, Const::PI);
        }

        inline void shuffle_nu() {
            if (orbit_type == OrbitType::ellipse) {
                Scalar M = Orbit::random_mean_anomaly(e, -Const::PI, Const::PI);
                Scalar E = Orbit::calc_eccentric_anomaly(M, e);
                nu = Orbit::calc_true_anomaly(E, e);
            } else {
                SPACEHUB_ABORT("Only elliptical orbit provides random anomaly method at this moment!");
            }
        }

        friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &obt) {
            SpaceH::display(os, obt.m1, obt.m2, obt.e, obt.p, obt.i, obt.Omega, obt.omega, obt.nu);
            return os;
        }
    };

    template<typename Particle, typename Scalar>
    void place_at_orbit(Particle &particle, OrbitArgs<Scalar> const &param) {
        using Vector = decltype(particle.pos);

        Scalar u = (param.m1 + param.m2) * Const::G;

        Scalar sin_nu = sin(param.nu);
        Scalar cos_nu = cos(param.nu);

        Scalar r = param.p / (1 + param.e * cos_nu);
        Scalar v = sqrt(u / param.p);

        particle.pos = r * Vector(cos_nu, sin_nu, 0);
        particle.vel = v * Vector(-sin_nu, param.e + cos_nu, 0);

        Orbit::euler_rotate(particle.pos, param.Omega, param.i, param.omega + Const::PI);
        Orbit::euler_rotate(particle.vel, param.Omega, param.i, param.omega + Const::PI);
    }

    template<typename Particle>
    struct Kepler {
        /* Typedef */
        using Scalar = typename Particle::Scalar;
        using Vector = typename Particle::Vector;
        using ParticleType = Particle;
    private:
        using Variant = std::variant<Scalar, RandomIndicator>;
    public:
        /* Typedef */

        Kepler() = delete;

        Kepler(OrbitArgs<Scalar> const &args) {
            if (args.m1 >= args.m2) {
                obj1_.mass = args.m1;
                obj2_.mass = args.m2;
            } else {
                obj1_.mass = args.m2;
                obj2_.mass = args.m1;
            }

            obj1_.pos = Vector(0.0, 0.0, 0.0);
            obj1_.vel = Vector(0.0, 0.0, 0.0);

            Scalar u = (args.m1 + args.m2) * Const::G;

            Scalar sin_nu = sin(args.nu);
            Scalar cos_nu = cos(args.nu);

            Scalar r = args.p / (1 + args.e * cos_nu);
            Scalar v = sqrt(u / args.p);

            obj2_.pos = r * Vector(cos_nu, sin_nu, 0);
            obj2_.vel = v * Vector(-sin_nu, args.e + cos_nu, 0);

            Orbit::euler_rotate(obj2_.pos, args.Omega, args.i, args.omega + Const::PI);
            Orbit::euler_rotate(obj2_.vel, args.Omega, args.i, args.omega + Const::PI);
        }

        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e, Variant i, Variant phi, Variant psi, Variant trueAnomaly)
                : Kepler(OrbitArgs<Scalar>(m1, m2, p, e, i, phi, psi, trueAnomaly)) {}

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


    template<typename Vector, typename Scalar>
    void coord_to_oribt_args(Scalar u, const Vector &dr, const Vector &dv, OrbitArgs<Scalar> &params) {
        Vector L = cross(dr, dv);
        Vector N = cross(Vector(0, 0, 1.0), L);
        Scalar r = norm(dr);
        Scalar n = norm(N);
        Scalar l = norm(L);
        Scalar rv = dot(dr, dv);
        Vector E = (dr * (norm2(dv) - u * re_norm(dr)) - dv * rv) / u;
        params.e = norm(E);

        if (!iseq(params.e, 1.0)) {
            Scalar a = -u * r / (r * norm2(dv) - 2.0 * u);
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
                params.nu = -SpaceH::sign(dot(N, dr)) * myacos(dot(peri / norm(peri), dr / r));
            } else {
                params.Omega = params.omega = 0;
                params.nu = SpaceH::sign(dr.y) * acos(dot(Vector(1.0, 0, 0), dr / r));
            }
        }
    }

    template<typename Vector, typename Scalar>
    Vector calc_eccentricity(Scalar u, Vector const &dr, Vector const &dv) {
        return (dr * (norm2(dv) - u * re_norm(dr)) - dv * dot(dr, dv)) / u;
    }

    template<typename Particle>
    auto calc_eccentricity(Particle const &p1, Particle const &p2) {
        return calc_eccentricity(Const::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
    }

    template<typename Vector, typename Scalar>
    Scalar calc_semi_major_axis(Scalar u, Vector const &dr, Vector const &dv) {
        Scalar r = norm(dr);
        return -u * r / (r * norm2(dv) - 2 * u);
    }

    template<typename Particle>
    auto calc_semi_major_axis(Particle const &p1, Particle const &p2) {
        return calc_semi_major_axis(Const::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
    }

    template<typename Vector, typename Scalar>
    auto calc_a_e(Scalar u, Vector const &dr, Vector const &dv) {
        Scalar r = norm(dr);
        Scalar v = norm(dv);
        Scalar v2 = v * v;
        Scalar vr = dot(dr, dv);
        Scalar vdfs = v2 - u / r;
        Scalar a = -u / (v2 - 2 * u / r);
        Scalar e = norm((dr * vdfs - dv * vr) / u);
        return std::make_tuple(a, e);
    }

    template<typename Particle>
    auto calc_a_e(Particle const &p1, Particle const &p2) {
        return calc_a_e(Const::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
    }
}
#endif
