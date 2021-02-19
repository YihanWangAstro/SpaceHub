/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file cross-section.hpp
 *
 * Header file.
 */
#pragma once

#include "../orbits/orbits.hpp"
#include "../orbits/particle-manip.hpp"
#include "../rand-generator.hpp"
#include "../vector/vector3.hpp"

/**
 * @namespace space::scattering
 * namespace for scattering
 */
namespace space::scattering {

    /**
     * @brief Calculate the critical velocity of the scattering between two clusters(can be single particle).
     *
     * @tparam Scalar Floating point like type.
     * @param[in] m1 Mass of the first cluster.
     * @param[in] m2 Mass of the second cluster.
     * @param[in] E1_inner The inner mechanical energy(mechanical energy in its own centre of mass frame) of the first
     * cluster.
     * @param[in] E2_inner The inner mechanical energy(mechanical energy in its own centre of mass frame) of the second
     * cluster.
     * @return The critical velocity.
     */
    template <typename Scalar>
    auto critical_vel(Scalar m1, Scalar m2, Scalar E1_inner, Scalar E2_inner) {
        auto m_rdc = m1 * m2 / (m1 + m2);
        return sqrt(-2 * (E1_inner + E2_inner) / m_rdc);
    }

    /**
     * @brief Calculate the critical velocity of the scattering between two clusters(can be single particle).
     * @tparam Cluster1 std::ranges(Container) with element type has public member `mass`(Scalar), `pos`(Vector) and
     * `vel`(Vector)./Type of single particle.
     * @tparam Cluster2 std::ranges(Container) with element type has public member `mass`(Scalar), `pos`(Vector) and
     * `vel`(Vector)./Type of single particle.
     * @param[in] stay_cluster The scattered cluster/scattered single particle.
     * @param[in] incident_cluster The incident cluster/incident single particle.
     * @return The critical velocity.
     */
    template <typename Cluster1, typename Cluster2>
    auto critical_vel(Cluster1 const& stay_cluster, Cluster2 const& incident_cluster) {
        auto M_stay = orbit::M_tot(stay_cluster);
        auto M_incident = orbit::M_tot(incident_cluster);
        auto E_inner1 = orbit::E_inner(stay_cluster);
        auto E_inner2 = orbit::E_inner(incident_cluster);
        return critical_vel(M_stay, M_incident, E_inner1, E_inner2);
    }

    /*template <typename Scalar>
    auto b_max(Scalar v_c, Scalar v_inf, Scalar a_max) {
      return a_max * (8 * v_c / v_inf + 3);
    }*/

    template <typename Scalar>
    auto b_max(Scalar m_tot, Scalar v_inf, Scalar rp_max) {
        return sqrt(rp_max * rp_max + 2 * m_tot * consts::G * rp_max / (v_inf * v_inf));
    }

    template <typename Cluster1, typename Cluster2, typename Scalar>
    auto b_max(Cluster1 const& stay_cluster, Cluster2 const& incident_cluster, Scalar v_inf, Scalar interact_factor) {
        auto const M_stay = orbit::M_tot(stay_cluster);
        auto const M_incident = orbit::M_tot(incident_cluster);

        auto const R1 = orbit::cluster_size(stay_cluster);
        auto const R2 = orbit::cluster_size(incident_cluster);

        auto const rp_max = orbit::tidal_radius(interact_factor, M_stay, M_incident, R1, R2);

        return b_max(M_stay + M_incident, v_inf, rp_max);
    }

    template <typename Scalar>
    auto b_max(Scalar M_stay, Scalar M_incident, Scalar R1, Scalar R2, Scalar v_inf, Scalar interact_factor) {
        auto const rp_max = orbit::tidal_radius(interact_factor, M_stay, M_incident, R1, R2);

        return b_max(M_stay + M_incident, v_inf, rp_max);
    }

    /**
     * @brief Randomly create an incident orbit that its infinity incident end is uniformly distributed in a circle area
     * with radius b_max.
     *
     * @tparam Scalar Floating point like type.
     * @param[in] m_stay The mass of the scattered object.
     * @param[in] m_incident The mass of the incident object.
     * @param[in] v_inf The relative velocity at infinity between scattering objects.
     * @param[in] b_max The max impact parameter.
     * @param[in] r The relative distance between the centre of mass between two objects. Used to generate the phase
     * (true anomaly of the orbit). The incident object will launch at this distance.
     * @return The incident hyperbolic orbit.
     */
    template <typename Scalar>
    auto incident_orbit(Scalar m_stay, Scalar m_incident, Scalar v_inf, Scalar b_max, Scalar r) {
        using Vector = Vec3<Scalar>;
        auto b = sqrt(random::Uniform(0, b_max * b_max));
        auto w = random::Uniform(0, 2 * consts::pi);
        return orbit::HyperOrbit(m_stay, m_incident, v_inf, b, w, 0.0, 0.0, r, orbit::Hyper::in);
    }

    /**
     * @brief Randomly create an incident orbit that its infinity incident end is uniformly distributed in a circle area
     * with radius b_max.
     * @tparam Cluster1 std::ranges(Container) with element type has public member `mass`(Scalar), `pos`(Vector) and
     * `vel`(Vector)./Type of single particle.
     * @tparam Cluster2 std::ranges(Container) with element type has public member `mass`(Scalar), `pos`(Vector) and
     * `vel`(Vector)./Type of single particle.
     * @tparam Scalar Floating point like type.
     * @param[in] stay_cluster The scattered cluster/scattered single particle.
     * @param[in] incident_cluster The incident cluster/incident single particle.
     * @param[in] v_inf The relative velocity at infinity between scattering objects.
     * @param[in] tidal_factor Tidal factor between two clusters. Used to generate the phase
     * (true anomaly of the orbit). The incident object will launch at the corresponding distance of this tidal factor.
     * @return The incident hyperbolic orbit.
     * TODO : b_max selection is implicit, need to be changed.
     * @deprecated
     */

    template <typename Cluster1, typename Cluster2, typename Scalar>
    auto incident_orbit(Cluster1 const& stay_cluster, Cluster2 const& incident_cluster, Scalar v_inf,
                        Scalar tidal_factor) {
        auto const M_stay = orbit::M_tot(stay_cluster);
        auto const M_incident = orbit::M_tot(incident_cluster);

        auto const R1 = orbit::cluster_size(stay_cluster);
        auto const R2 = orbit::cluster_size(incident_cluster);

        Scalar interact_factor = 0.02;
        auto const closest_approach_max = orbit::tidal_radius(interact_factor, M_stay, M_incident, R1, R2);

        auto const b_upper = b_max(M_stay + M_incident, v_inf, closest_approach_max);

        auto const r_start = orbit::tidal_radius(tidal_factor, M_stay, M_incident, R1, R2);
        return incident_orbit(M_stay, M_incident, v_inf, b_upper, r_start);
    }

    template <typename Cluster1, typename Cluster2, typename Scalar>
    auto incident_orbit(Cluster1 const& stay_cluster, Cluster2 const& incident_cluster, Scalar v_inf, Scalar b_max,
                        Scalar tidal_factor) {
        auto const M_stay = orbit::M_tot(stay_cluster);
        auto const M_incident = orbit::M_tot(incident_cluster);
        auto const R1 = orbit::cluster_size(stay_cluster);
        auto const R2 = orbit::cluster_size(incident_cluster);

        auto const r_start = orbit::tidal_radius(tidal_factor, M_stay, M_incident, R1, R2);

        return incident_orbit(M_stay, M_incident, v_inf, b_max, r_start);
    }

    template <typename Scalar>
    inline auto hard_radius(Scalar m1, Scalar m2, Scalar m_evn, Scalar sigma) {
        return consts::G * m1 * m2 / (m_evn * sigma * sigma);
    }

}  // namespace space::scattering
