#ifndef ORBITS_H
#define ORBITS_H

#include "src/type_class.h"
#include "src/own_math.h"
#include "../macros.h"
#include <math.h>

namespace SpaceH {
    namespace Orbits {
        /**
         *
         * @tparam Scalar
         * @param x
         * @return
         */
        template<typename Scalar>
        inline Scalar myacos(Scalar x) {
            return acos(SpaceH::min(SpaceH::max(-1.0, x), 1.0));
        }

        /**
         *
         * @tparam Scalar
         */
        template<typename Scalar>
        struct Param {
            Scalar e;
            Scalar a;
            Scalar p;
            Scalar i;
            Scalar Omega;
            Scalar omega;
            Scalar nu;//true anomaly
            friend std::ostream &operator<<(std::ostream &os, const Param &P) {
                os << P.e << ' ' << P.a << ' ' << P.p << ' ' << P.i << ' '
                   << P.Omega << ' ' << P.omega << ' ' << P.nu;
                return os;
            }
        };

        /**
         *
         * @tparam Vector
         * @tparam Scalar
         * @param m1
         * @param m2
         * @param pos1
         * @param pos2
         * @param vel1
         * @param vel2
         * @return
         */
        template<typename Vector, typename Scalar>
        inline Vector
        calcuEccentricity(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                          const Vector &vel2) {
            Vector dr = pos1 - pos2;
            Vector dv = vel1 - vel2;
            Scalar u = Const::G * (m1 + m2);
            return (dr * (dv.norm2() - u * dr.reNorm()) - dv * dot(dr, dv)) / u;
        }

        /**
         *
         * @tparam Vector
         * @tparam Scalar
         * @param m1
         * @param m2
         * @param pos1
         * @param pos2
         * @param vel1
         * @param vel2
         * @return
         */
        template<typename Vector, typename Scalar>
        inline Scalar
        calcuSemiMajorAxis(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                           const Vector &vel2) {
            Vector dr = pos1 - pos2;
            Vector dv = vel1 - vel2;
            Scalar u = Const::G * (m1 + m2);
            Scalar r = dr.norm();
            return -u * r / (r * dv.norm2() - 2 * u);
        }

        /**
         *
         * @tparam Vector
         * @tparam Scalar
         * @param m1
         * @param m2
         * @param pos1
         * @param pos2
         * @param vel1
         * @param vel2
         * @param params
         */
        template<typename Vector, typename Scalar>
        inline void toOribtParameters(Scalar m1, Scalar m2, const Vector &pos1, const Vector &pos2, const Vector &vel1,
                                      const Vector &vel2, Param<Scalar> &params) {
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

            if (fabs(params.e - 1.0) > SpaceH::epsilon<Scalar>::value) {
                params.a = -u * r / (r * dv.norm2() - 2.0 * u);
                params.p = params.a * (1.0 - params.e * params.e);
            } else {
                params.a = 1.0 / 0;
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

        /**
         *
         * @tparam Vector
         * @tparam Scalar
         * @param v
         * @param phi
         * @param theta
         * @param psi
         */
        template<typename Vector, typename Scalar>
        void eulerRotate(Vector &v, const Scalar phi, const Scalar theta, const Scalar psi) {
            Vector tmp(0.0, 0.0, 0.0);

            Scalar cos_phi   = cos(phi);
            Scalar sin_phi   = sin(phi);
            Scalar cos_psi   = cos(psi);
            Scalar sin_psi   = sin(psi);
            Scalar cos_theta = cos(theta);
            Scalar sin_theta = sin(theta);

            tmp.x =   v.x * (cos_phi * cos_psi - sin_phi * cos_theta * sin_psi)
                    - v.y * (cos_phi * sin_psi + sin_phi * cos_theta * cos_psi)
                    + v.z * (sin_phi * sin_theta);

            tmp.y =   v.x * (sin_phi * cos_psi + cos_phi * cos_theta * sin_psi)
                    - v.y * (sin_phi * sin_psi - cos_phi * cos_theta * cos_psi)
                    - v.z * (cos_phi * sin_theta);

            tmp.z = v.x * sin_theta * sin_psi + v.y * sin_theta * cos_psi + v.z * cos_theta;

            v = tmp;
        };

        /**
         *
         * @tparam Scalar
         * @param e
         * @param Mmin
         * @param Mmax
         * @return
         */
        template<typename Scalar>
        Scalar getRandomMeanAnomaly(Scalar e, Scalar Mmin, Scalar Mmax) {
            SpaceH::Random<T> uniform;
            if (e >= 0)
                return uniform() * (Mmax - Mmin) + Mmin;
            else {
                SpaceH::ERR_MSG("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }

        /**
         *
         * @tparam Scalar
         * @param E
         * @param e
         * @return
         */
        template<typename Scalar>
        Scalar calcuTrueAnomaly(Scalar E, Scalar e) {
            if (0 <= e && e < 1)
                return 2 * atan2(sqrt(1 + e) * sin(E * 0.5), sqrt(1 - e) * cos(0.5 * E));
            else if (e > 1)
                return 2 * atan2(sqrt(1 + e) * sinh(E * 0.5), sqrt(e - 1) * cosh(0.5 * E));
            else if (fabs(e - 1) < SpaceH::epsilon<Scalar>::value)
                return 2 * atan(0.5 * E);
            else {
                SpaceH::ERR_MSG("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }

        /**
         *
         * @tparam Scalar
         * @param M
         * @param e
         * @return
         */
        template<typename Scalar>
        Scalar calcuEccentricAnomaly(Scalar M, Scalar e) {
            if (0 <= e && e < 1)//newton iteration may encounter stationary point
                return SpaceH::root_dichotom(
                        [&](Scalar x) -> Scalar { return x - e * sin(x) - M; });//find this function in ownMath.h
            else if (e > 1)
                return SpaceH::root_dichotom([&](Scalar x) -> Scalar { return e * sinh(x) - x - M; });
            else if (fabs(e - 1) < SpaceH::epsilon<Scalar>::value)
                return SpaceH::root_dichotom([&](Scalar x) -> Scalar { return x + x * x * x / 3 - M; });
            else {
                SpaceH::ERR_MSG("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }
    }

    /**
     *
     * @tparam Dtype
     */
    template<typename Dtype>
    struct Particle {
        /* Typedef */
        using type   = SpaceH::ProtoType<Dtype, SpaceH::DYNAMICAL>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */

        Vector pos;
        Vector vel;
        Scalar mass;
        Scalar radius;

        friend std::ostream &operator<<(std::ostream &os, const Particle &P) {
            os << P.pos / Unit::AU << ' ' << P.vel / Unit::KMS << ' ' << P.mass / Unit::M_SOLAR << ' '
               << P.radius / Unit::AU;
            return os;
        }
    };

    /**
     *
     * @tparam Dtype
     */
    template<typename Dtype>
    struct Kepler {
        /* Typedef */
        using type   = SpaceH::ProtoType<Dtype, SpaceH::DYNAMICAL>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;

        /* Typedef */

        inline const Scalar a() const { return param.a; }

        inline const Scalar b() const { return b_; }

        inline const Scalar p() const { return param.p; }//semi-latus rectum

        inline const Scalar e() const { return param.e; }

        inline const Scalar nu() const { return param.nu; }

        inline const Scalar Omega() const { return param.Omega; }

        inline const Scalar i() const { return param.i; }

        inline const Scalar omega() const { return param.omega; }

        Kepler() = delete;

        /**
         *
         * @param m1
         * @param m2
         * @param p
         * @param e
         */
        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e) {
            SpaceH::Random<T> uniform;
            calcuOrbitalParameter(m1, m2, p, e);
            param.Omega = uniform() * 2 * Const::PI - Const::PI;
            param.i     = acos(SpaceH::uniform() * 2 - 1);
            param.omega = uniform() * 2 * Const::PI - Const::PI;
            Scalar M    = Orbits::getRandomMeanAnomaly(e, -Unit::HUBBLETIME / T_, Unit::HUBBLETIME / T_);
            Scalar E    = Orbits::calcuEccentricAnomaly(M, param.e);
            param.nu    = Orbits::calcuTrueAnomaly(E, param.e);
            createOrbit(m1, m2);
        }

        /**
         *
         * @param m1
         * @param m2
         * @param p
         * @param e
         * @param phi
         * @param inclination
         * @param psi
         * @param trueAnomaly
         */
        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e, Scalar phi, Scalar inclination, Scalar psi,
               Scalar trueAnomaly = Const::PI) {
            calcuOrbitalParameter(m1, m2, p, e);
            param.Omega = phi;
            param.i     = inclination;
            param.omega = psi;
            param.nu    = trueAnomaly;
            createOrbit(m1, m2);
        }

        /**
         *
         * @param newCMPos
         * @param newCMVel
         */
        void moveOrbitTo(const Vector &newCMPos, const Vector &newCMVel) {
            moveToCentreMassCoord();
            P1.pos += newCMPos, P1.vel += newCMVel;
            P2.pos += newCMPos, P2.vel += newCMVel;
        }

        /**
         *
         * @param P
         */
        void moveOrbitTo(const Particle<Scalar> &P) {
            moveToCentreMassCoord();
            P1.pos += P.pos, P1.vel += P.vel;
            P2.pos += P.pos, P2.vel += P.vel;
        }

        /**
         *
         */
        void moveToCentreMassCoord() {
            Vector CMP = getCentreMassPos();
            Vector CMV = getCentreMassVel();
            P1.pos -= CMP, P1.vel -= CMV;
            P2.pos -= CMP, P2.vel -= CMV;
        }

        /**
         *
         */
        void moveToPrimaryCentre() {
            if (P1.mass > P2.mass) {
                P2.pos -= P1.pos, P1.pos.setZero();
                P2.vel -= P1.vel, P1.vel.setZero();
            } else {
                P1.pos -= P2.pos, P2.pos.setZero();
                P1.vel -= P2.vel, P2.vel.setZero();
            }
        }

        /**
         *
         * @return
         */
        inline Vector getCentreMassPos() {
            Scalar totalMass = P1.mass + P2.mass;
            return Vector((P1.mass * P1.pos + P2.mass * P2.pos) / totalMass);
        }

        /**
         *
         * @return
         */
        inline Vector getCentreMassVel() {
            Scalar totalMass = P1.mass + P2.mass;
            return Vector((P1.mass * P1.vel + P2.mass * P2.vel) / totalMass);
        }

        /**
         *
         * @return
         */
        inline const Particle<Dtype> &primary() const {
            if (P1.mass >= P2.mass)
                return P1;
            else
                return P2;
        }

        /**
         *
         * @return
         */
        inline const Particle<Dtype> &secondary() const {
            if (P1.mass < P2.mass)
                return P1;
            else
                return P2;
        }

        /**
         *
         */
        void randomPhi() {
            SpaceH::Random<T> uniform;
            param.Omega = uniform() * 2 * Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }

        /**
         *
         */
        void randomTheta() {
            SpaceH::Random<T> uniform;
            param.i = acos(uniform() * 2 - 1);
            createOrbit(P1.mass, P2.mass);
        }

        /**
         *
         */
        void randomPsi() {
            SpaceH::Random<T> uniform;
            param.omega = uniform() * 2 * Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }

        /**
         *
         */
        void randomAngles() {
            SpaceH::Random<T> uniform;
            param.Omega = uniform() * 2 * Const::PI - Const::PI;
            param.i = acos(uniform() * 2 - 1);
            param.omega = uniform() * 2 * Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }

        /**
         *
         * @param pastTime
         * @param futureTime
         */
        void randomPhase(Scalar pastTime = -Unit::HUBBLETIME, Scalar futureTime = Unit::HUBBLETIME) {
            Scalar M = 0;
            if (param.e < 1 && futureTime - pastTime > T_)
                M = Orbits::getRandomMeanAnomaly(param.e, -SpaceH::Const::PI, SpaceH::Const::PI);
            else
                M = Orbits::getRandomMeanAnomaly(param.e, pastTime / T_, futureTime / T_);
            Scalar E = Orbits::calcuEccentricAnomaly(M, param.e);
            param.nu = Orbits::calcuTrueAnomaly(E, param.e);
            createOrbit(P1.mass, P2.mass);
        }

    private:
        /**
         *
         * @param m1
         * @param m2
         * @param p
         * @param e
         */
        void calcuOrbitalParameter(Scalar m1, Scalar m2, Scalar p, Scalar e) {
            checkParameter(p, e);
            P1.mass = m1 * Unit::M_SOLAR;
            P2.mass = m2 * Unit::M_SOLAR;
            param.p = p * Unit::AU;
            param.e = e;
            if (fabs(param.e - 1) > SpaceH::epsilon<Scalar>::value) {
                param.a = param.p / (1 - param.e * param.e);
                b_ = fabs(param.a) * sqrt(fabs(1 - param.e * param.e));
                T_ = sqrt(fabs(param.a * param.a * param.a / Const::G / (P1.mass + P2.mass)));
            } else {
                param.a = 1.0 / 0.0;
                b_ = 0.5 * param.p;
                T_ = sqrt(fabs(0.25 * param.p * param.p * param.p / Const::G / (P1.mass + P2.mass)));
            }
        }

        /**
         *
         * @param p
         * @param e
         */
        void checkParameter(Scalar p, Scalar e) {
            if (p < 0)
                SpaceH::ERR_MSG("semi-latus rectum cannot be negative", __FILE__, __LINE__);
            if (e < 0)
                SpaceH::ERR_MSG("Eccentrcity cannot be negative!", __FILE__, __LINE__);
        }

        /**
         *
         * @param m1
         * @param m2
         */
        void createOrbit(Scalar m1, Scalar m2) {
            Scalar u = (m1 + m2) * Const::G;

            Scalar r = param.p / (1 + param.e * cos(param.nu));
            P2.pos = r * Vector(cos(param.nu), sin(param.nu), 0);

            Scalar v = sqrt(u / param.p);
            P2.vel = v * Vector(-sin(param.nu), param.e + cos(param.nu), 0);

            P1.pos.setZero(), P1.vel.setZero();

            P1.mass = m1 * Unit::M_SOLAR, P1.radius = 0;
            P2.mass = m2 * Unit::M_SOLAR, P2.radius = 0;

            Orbits::eulerRotate(P2.pos, param.Omega, param.i, param.omega + Const::PI);
            Orbits::eulerRotate(P2.vel, param.Omega, param.i, param.omega + Const::PI);
            moveToCentreMassCoord();
        }

    private:
        Particle<Dtype> P1;
        Particle<Dtype> P2;
        Orbits::Param<Scalar> param;
        Scalar b_;
        Scalar T_;
    };
}
#endif
