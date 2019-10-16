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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file symplectic-integrator.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_SYMPLECTIC_INTEGRATOR_HPP
#define SPACEHUB_SYMPLECTIC_INTEGRATOR_HPP

#include "../../particle-system/particle-system.hpp"
/**
 * @namespace space::integrator
 * Documentation for space
 */
namespace space::integrator {

  /**
   *
   * @tparam Derived
   */
  template<typename Derived>
  class SymIntegrator {
  public:
    static constexpr size_t order{Derived::order};

    template<typename T>
    void integrate(T &particle_system, typename T::Scalar step_size) {
      static_assert(particle_system::is_particle_system_v<T>, "Passing non paritcle-system-type!");
      static_cast<Derived *>(this)->impl_integrate(particle_system, step_size);
    }

  private:
    SymIntegrator() = default;

    friend Derived;
  };

  /**
   *
   */
  class symplectic2nd : public SymIntegrator<symplectic2nd> {
  public:
    static constexpr size_t order{2};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &system, typename ParticleSys::Scalar step_size) {
      system.drift(0.5 * step_size);
      system.kick(step_size);
      system.drift(0.5 * step_size);
    }
  };

  /**
   *
   */
  class symplectic4th : public SymIntegrator<symplectic4th> {
  public:
    static constexpr size_t order{4};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &system, typename ParticleSys::Scalar step_size) {
      system.drift(6.7560359597983000E-1 * step_size);
      system.kick(1.3512071919596600E0 * step_size);
      system.drift(-1.7560359597983000E-1 * step_size);
      system.kick(-1.7024143839193200E0 * step_size);
      system.drift(-1.7560359597983000E-1 * step_size);
      system.kick(1.3512071919596600E0 * step_size);
      system.drift(6.7560359597983000E-1 * step_size);
    }
  };

  /**
   *
   */
  class symplectic6th : public SymIntegrator<symplectic6th> {
  public:
    static constexpr size_t order{6};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &system, typename ParticleSys::Scalar step_size) {
      /*unroll loop manually*/
      system.drift(3.9225680523877998E-1 * step_size);
      system.kick(7.8451361047755996E-1 * step_size);
      system.drift(5.1004341191845848E-1 * step_size);
      system.kick(2.3557321335935699E-1 * step_size);
      system.drift(-4.7105338540975655E-1 * step_size);
      system.kick(-1.1776799841788701E0 * step_size);
      system.drift(6.8753168252518093E-2 * step_size);
      system.kick(1.3151863206839063E0 * step_size);
      system.drift(6.8753168252518093E-2 * step_size);
      system.kick(-1.1776799841788701E0 * step_size);
      system.drift(-4.7105338540975655E-1 * step_size);
      system.kick(2.3557321335935699E-1 * step_size);
      system.drift(5.1004341191845848E-1 * step_size);
      system.kick(7.8451361047755996E-1 * step_size);
      system.drift(3.9225680523877998E-1 * step_size);
    }
  };

  /**
   *
   */
  class symplectic8th : public SymIntegrator<symplectic8th> {
  public:
    static constexpr size_t order{8};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &system, typename ParticleSys::Scalar step_size) {
      /*unroll loop manually*/
      system.drift(5.21213104349955048E-1 * step_size);
      system.kick(1.04242620869991010E0 * step_size);
      system.drift(1.43131625920352512E0 * step_size);
      system.kick(1.82020630970713992E0 * step_size);
      system.drift(9.88973118915378424E-1 * step_size);
      system.kick(1.57739928123617007E-1 * step_size);
      system.drift(1.29888362714548355E0 * step_size);
      system.kick(2.44002732616735019E0 * step_size);
      system.drift(1.21642871598513458E0 * step_size);
      system.kick(-7.16989419708119989E-3 * step_size);
      system.drift(-1.22708085895116059E0 * step_size);
      system.kick(-2.44699182370524015E0 * step_size);
      system.drift(-2.03140778260310517E0 * step_size);
      system.kick(-1.61582374150096997E0 * step_size);
      system.drift(-1.69832618404521085E0 * step_size);
      system.kick(-1.78082862658945151E0 * step_size);
      system.drift(-1.69832618404521085E0 * step_size);
      system.kick(-1.61582374150096997E0 * step_size);
      system.drift(-2.03140778260310517E0 * step_size);
      system.kick(-2.44699182370524015E0 * step_size);
      system.drift(-1.22708085895116059E0 * step_size);
      system.kick(-7.16989419708119989E-3 * step_size);
      system.drift(1.21642871598513458E0 * step_size);
      system.kick(2.44002732616735019E0 * step_size);
      system.drift(1.29888362714548355E0 * step_size);
      system.kick(1.57739928123617007E-1 * step_size);
      system.drift(9.88973118915378424E-1 * step_size);
      system.kick(1.82020630970713992E0 * step_size);
      system.drift(1.43131625920352512E0 * step_size);
      system.kick(1.04242620869991010E0 * step_size);
      system.drift(5.21213104349955048E-1 * step_size);
    }
  };

  /**
   *
   */
  class symplectic10th : public SymIntegrator<symplectic10th> {
  public:
    static constexpr size_t order{10};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &system, typename ParticleSys::Scalar step_size) {
      /*unroll loop manually*/
      system.drift(3.0610967201933609e-01 * step_size);
      system.kick(6.1221934403867218e-01 * step_size);
      system.drift(-9.4012698954724694e-02 * step_size);
      system.kick(-8.0024474194812156e-01 * step_size);
      system.drift(-6.6002635995076209e-01 * step_size);
      system.kick(-5.1980797795340250e-01 * step_size);
      system.drift(-1.5240397828727220e-01 * step_size);
      system.kick(2.1500002137885812e-01 * step_size);
      system.drift(-1.1750569210727700e-01 * step_size);
      system.kick(-4.5001140559341213e-01 * step_size);
      system.drift(2.2250778443570857e-01 * step_size);
      system.kick(8.9502697446482926e-01 * step_size);
      system.drift(5.1288848042847668e-01 * step_size);
      system.kick(1.3074998639212410e-01 * step_size);
      system.drift(3.3095796002497074e-01 * step_size);
      system.kick(5.3116593365781739e-01 * step_size);
      system.drift(-6.0050191119721985e-02 * step_size);
      system.kick(-6.5126631589726136e-01 * step_size);
      system.drift(-7.6956706144236287e-01 * step_size);
      system.kick(-8.8786780698746448e-01 * step_size);
      system.drift(-7.6872229417056015e-02 * step_size);
      system.kick(7.3412334815335245e-01 * step_size);
      system.drift(4.2477286784491525e-01 * step_size);
      system.kick(1.1542238753647800e-01 * step_size);
      system.drift(4.3160892192959932e-01 * step_size);
      system.kick(7.4779545632272060e-01 * step_size);
      system.drift(5.5434862753225678e-02 * step_size);
      system.kick(-6.3692573081626924e-01 * step_size);
      system.drift(-1.9288621063874828e-01 * step_size);
      system.kick(2.5115330953877268e-01 * step_size);
      system.drift(3.3904387248169282e-01 * step_size);
      system.kick(4.2693443542461296e-01 * step_size);
      system.drift(3.3904387248169282e-01 * step_size);
      system.kick(2.5115330953877268e-01 * step_size);
      system.drift(-1.9288621063874828e-01 * step_size);
      system.kick(-6.3692573081626924e-01 * step_size);
      system.drift(5.5434862753225678e-02 * step_size);
      system.kick(7.4779545632272060e-01 * step_size);
      system.drift(4.3160892192959932e-01 * step_size);
      system.kick(1.1542238753647800e-01 * step_size);
      system.drift(4.2477286784491525e-01 * step_size);
      system.kick(7.3412334815335245e-01 * step_size);
      system.drift(-7.6872229417056015e-02 * step_size);
      system.kick(-8.8786780698746448e-01 * step_size);
      system.drift(-7.6956706144236287e-01 * step_size);
      system.kick(-6.5126631589726136e-01 * step_size);
      system.drift(-6.0050191119721985e-02 * step_size);
      system.kick(5.3116593365781739e-01 * step_size);
      system.drift(3.3095796002497074e-01 * step_size);
      system.kick(1.3074998639212410e-01 * step_size);
      system.drift(5.1288848042847668e-01 * step_size);
      system.kick(8.9502697446482926e-01 * step_size);
      system.drift(2.2250778443570857e-01 * step_size);
      system.kick(-4.5001140559341213e-01 * step_size);
      system.drift(-1.1750569210727700e-01 * step_size);
      system.kick(2.1500002137885812e-01 * step_size);
      system.drift(-1.5240397828727220e-01 * step_size);
      system.kick(-5.1980797795340250e-01 * step_size);
      system.drift(-6.6002635995076209e-01 * step_size);
      system.kick(-8.0024474194812156e-01 * step_size);
      system.drift(-9.4012698954724694e-02 * step_size);
      system.kick(6.1221934403867218e-01 * step_size);
      system.drift(3.0610967201933609e-01 * step_size);
    }
  };

  template<typename T>
  constexpr bool is_sym_integrator_v = std::is_base_of_v<SymIntegrator<T>, T>;
}
#endif //SPACEHUB_SYMPLECTIC_INTEGRATOR_HPP
