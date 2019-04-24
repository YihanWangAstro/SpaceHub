#ifndef ORBITS_H
#define ORBITS_H

#include "../rand-generator.tpp"
#include "../own-math.h"
#include "../macros.h"
#include "../core-computation.tpp"
#include <math.h>
#include <variant>

namespace space::orbit {

    template<typename Scalar>
    inline Scalar myacos(Scalar x) {
        return acos(space::min(space::max(-1.0, x), 1.0));
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
    Scalar calc_true_anomaly(Scalar E, Scalar e) {
        if (0 <= e && e < 1)
            return 2 * atan2(sqrt(1 + e) * sin(E * 0.5), sqrt(1 - e) * cos(0.5 * E));
        else if (e > 1)
            return 2 * atan2(sqrt(1 + e) * sinh(E * 0.5), sqrt(e - 1) * cosh(0.5 * E));
        else if (iseq(e, 1.0))
            return 2 * atan(0.5 * E);
        else {
            spacehub_abort("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }

    template<typename Scalar>
    Scalar calc_eccentric_anomaly(Scalar M, Scalar e) {
        if (0 <= e && e < 1)//newton iteration may encounter stationary point
            return space::root_dichotom(
                    [&](Scalar x) -> Scalar { return x - e * sin(x) - M; });//find this function in ownMath.h
        else if (e > 1)
            return space::root_dichotom([&](Scalar x) -> Scalar { return e * sinh(x) - x - M; });
        else if (fabs(e - 1) < space::epsilon<Scalar>::value)
            return space::root_dichotom([&](Scalar x) -> Scalar { return x + x * x * x / 3 - M; });
        else {
            spacehub_abort("Eccentrcity cannot be negative, Nan or inf!");
            return 0;
        }
    }

    enum class OrbitType {
        ellipse, parabola, hyperbola, none
    };

    template<typename T>
    constexpr OrbitType classify_orbit(T eccentricity) {
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
        Scalar u;//gravitational parameter
        Scalar e;//eccentricity
        Scalar p;//semi-latus rectum
        Scalar i;//inclination
        Scalar Omega;//longitude of the ascending node
        Scalar omega;//argument of periapsis
        Scalar nu;//true anomaly
        OrbitType orbit_type;

        OrbitArgs() = delete;

        OrbitArgs(Scalar tot_mass, Scalar _p_, Scalar _e_, Variant tilt, Variant LoAN, Variant AoP,
                  Variant true_anomaly) {
            if (_p_ < 0) spacehub_abort("Semi-latus rectum cannot be negative");

            orbit_type = classify_orbit(_e_);

            if (orbit_type == OrbitType::none) spacehub_abort("Eccentrcity cannot be negative or NaN!");

            u = tot_mass * consts::G;
            p = _p_;
            e = _e_;

            if (std::holds_alternative<Scalar>(tilt)) {
                i = std::get<Scalar>(tilt);
            } else {
                shuffle_i();
            }

            if (std::holds_alternative<Scalar>(LoAN)) {
                Omega = std::get<Scalar>(LoAN);
            } else {
                shuffle_Omega();
            }

            if (std::holds_alternative<Scalar>(AoP)) {
                omega = std::get<Scalar>(AoP);
            } else {
                shuffle_omega();
            }

            if (std::holds_alternative<Scalar>(true_anomaly)) {
                nu = std::get<Scalar>(true_anomaly);
            } else {
                shuffle_nu();
            }
        }

        inline void shuffle_i() {
            i = acos(randomGen::Uniform<Scalar>::get(-1, 1));
        }

        inline void shuffle_Omega() {
            Omega = randomGen::Uniform<Scalar>::get(-consts::pi, consts::pi);
        }

        inline void shuffle_omega() {
            omega = randomGen::Uniform<Scalar>::get(-consts::pi, consts::pi);
        }

        inline void shuffle_nu() {
            if (orbit_type == OrbitType::ellipse) {
                Scalar M = randomGen::Uniform<Scalar>::get(-consts::pi, consts::pi);
                Scalar E = orbit::calc_eccentric_anomaly(M, e);
                nu = orbit::calc_true_anomaly(E, e);
            } else {
                spacehub_abort("Only elliptical orbit provides random anomaly method at this moment!");
            }
        }

        friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &obt) {
            space::display(os, obt.u, obt.e, obt.p, obt.i, obt.Omega, obt.omega, obt.nu);
            return os;
        }
    };

    using Kepler = OrbitArgs<double>;

    template<typename Vector, typename Scalar>
    void coord_to_oribt_args(Scalar u, const Vector &dr, const Vector &dv, OrbitArgs<Scalar> &args) {
        Vector L = cross(dr, dv);
        Vector N = cross(Vector(0, 0, 1.0), L);
        Scalar r = norm(dr);
        Scalar n = norm(N);
        Scalar l = norm(L);
        Scalar rv = dot(dr, dv);
        Vector E = (dr * (norm2(dv) - u * re_norm(dr)) - dv * rv) / u;

        args.u = u;
        args.e = norm(E);
        args.orbit_type = classify_orbit(args.e);

        if (args.orbit_type == OrbitType::parabola) {
            Scalar a = -u * r / (r * norm2(dv) - 2.0 * u);
            args.p = semi_latus_rectum(a, args.e);
        } else {
            args.p = l * l / u;
        }

        args.i = acos(L.z / l);

        if (args.e != 0) {
            args.nu = space::sign(rv) * myacos(dot(E / args.e, dr / r));

            if (n != 0) {
                args.Omega = space::sign(N.y) * myacos(N.x / n);
                args.omega = space::sign(E.z) * myacos(dot(E / args.e, N / n));
            } else {
                args.omega = -space::sign(E.y) * myacos(-E.x / args.e);
                args.Omega = args.omega;
            }
        } else {
            if (n != 0) {
                args.Omega = space::sign(N.y) * myacos(N.x / n);
                args.omega = 0;
                Vector peri = cross(L, N);
                args.nu = -space::sign(dot(N, dr)) * myacos(dot(peri / norm(peri), dr / r));
            } else {
                args.Omega = args.omega = 0;
                args.nu = space::sign(dr.y) * acos(dot(Vector(1.0, 0, 0), dr / r));
            }
        }
    }

    template<typename Vector, typename Scalar>
    void oribt_args_to_coord(OrbitArgs<Scalar> const &args, Vector &pos, Vector &vel) {
        Scalar u = args.u;

        Scalar sin_nu = sin(args.nu);
        Scalar cos_nu = cos(args.nu);

        Scalar r = args.p / (1 + args.e * cos_nu);
        Scalar v = sqrt(u / args.p);

        pos = r * Vector(cos_nu, sin_nu, 0);
        vel = v * Vector(-sin_nu, args.e + cos_nu, 0);

        orbit::euler_rotate(pos, args.Omega, args.i, args.omega + consts::pi);
        orbit::euler_rotate(vel, args.Omega, args.i, args.omega + consts::pi);
    }

    template<typename Particle, typename ...Args>
    void move_to_com_coord(Particle &ptc, Args &...ptcs) {
        if constexpr (sizeof...(Args) != 0) {
            static_assert(calc::all(std::is_same_v<Args, Particle>...), "Wrong particle type!");
            auto tot_mass = (ptcs.mass + ... + ptc.mass);
            auto cm_pos = ((ptcs.mass * ptcs.pos) + ... + (ptc.mass * ptc.pos)) / tot_mass;
            auto cm_vel = ((ptcs.mass * ptcs.vel) + ... + (ptc.mass * ptc.vel)) / tot_mass;
            ((ptc.pos -= cm_pos), ..., (ptcs.pos -= cm_pos));
            ((ptc.vel -= cm_vel), ..., (ptcs.vel -= cm_vel));
        } else if constexpr (is_container_v<Particle>) {
            using sParticle = typename Particle::value_type;
            using Scalar = typename sParticle::Scalar;
            using Vector = typename sParticle::Vector;

            Scalar tot_mass{0};
            Vector cm_pos{0};
            Vector cm_vel{0};

            for (auto const &p : ptc) {
                tot_mass += p.mass;
                cm_pos += p.mass * p.pos;
                cm_vel += p.mass * p.vel;
            }

            cm_pos /= tot_mass;
            cm_vel /= tot_mass;

            for (auto &p : ptc) {
                p.pos -= cm_pos;
                p.vel -= cm_vel;
            }
        }
    }

    template<typename Vector, typename Particle, typename ...Args>
    void move_particles_to(Vector const &cm_pos, Vector const &cm_vel, Particle &ptc, Args &...ptcs) {
        if constexpr (sizeof...(Args) != 0) {
            static_assert(calc::all(std::is_same_v<Args, Particle>...), "Wrong particle type!");
            move_to_com_coord(ptc, ptcs...);
            ((ptc.pos += cm_pos), ..., (ptcs.pos += cm_pos));
            ((ptc.vel += cm_vel), ..., (ptcs.vel += cm_vel));
        } else if constexpr (is_container_v<Particle>) {
            move_to_com_coord(ptc);
            for (auto &p : ptc) {
                p.pos += cm_pos;
                p.vel += cm_vel;
            }
        } else {
            ptc.pos = cm_pos;
            ptc.vel = cm_vel;
        }
    }

    template<typename Particle, typename ...Args>
    void move_particles_to(Kepler const &args, Particle &ptc, Args &...ptcs) {
        using Vector = typename Particle::Vector;
        Vector cm_pos, cm_vel;
        oribt_args_to_coord(args, cm_pos, cm_vel);
        move_particles_to(cm_pos, cm_vel, ptc, ptcs...);
    }

    template<typename Scalar>
    inline constexpr auto semi_latus_rectum(Scalar a, Scalar e) {
        return a * (1 - e * e);
    }

    template<typename Vector, typename Scalar>
    Vector calc_eccentricity(Scalar u, Vector const &dr, Vector const &dv) {
        return (dr * (norm2(dv) - u * re_norm(dr)) - dv * dot(dr, dv)) / u;
    }

    template<typename Particle>
    auto calc_eccentricity(Particle const &p1, Particle const &p2) {
        return calc_eccentricity(consts::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
    }

    template<typename Vector, typename Scalar>
    Scalar calc_semi_major_axis(Scalar u, Vector const &dr, Vector const &dv) {
        Scalar r = norm(dr);
        return -u * r / (r * norm2(dv) - 2 * u);
    }

    template<typename Particle>
    auto calc_semi_major_axis(Particle const &p1, Particle const &p2) {
        return calc_semi_major_axis(consts::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
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
        return calc_a_e(consts::G * (p1.mass + p2.mass), p1.pos - p2.pos, p1.vel - p2.vel);
    }
}
#endif
