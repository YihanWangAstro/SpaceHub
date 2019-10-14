//
// Created by root on 10/7/19.
//
#include "dev-tools.hpp"

#ifndef SPACEHUB_ACCELERATIONS_H
#define SPACEHUB_ACCELERATIONS_H
namespace space {
/*---------------------------------------------------------------------------*\
    Class Accelerations Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Interactions
 * @tparam Coord
 */
  template<typename Interactions, typename Coord>
  class Accelerations {
  public:
    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(Accelerations, default, default, default, default, default);

    explicit Accelerations(size_t size) : acc_{size}, newtonian_acc_{size}, tot_vel_indep_acc_{size} {
      if constexpr (Interactions::ext_vel_indep) {
        ext_vel_indep_acc_.resize(size);
      }
      if constexpr (Interactions::ext_vel_dep) {
        ext_vel_dep_acc_.resize(size);
      }
    };
    
    //Public methods
    SPACEHUB_STD_ACCESSOR(auto, acc, acc_);

    SPACEHUB_STD_ACCESSOR(auto, newtonian_acc, newtonian_acc_);

    SPACEHUB_STD_ACCESSOR(auto, tot_vel_indep_acc, tot_vel_indep_acc_);

    SPACEHUB_STD_ACCESSOR(auto, ext_vel_indep_acc, ext_vel_indep_acc_);

    SPACEHUB_STD_ACCESSOR(auto, ext_vel_dep_acc, ext_vel_dep_acc_);

  private:
    Coord acc_;

    Coord newtonian_acc_;

    Coord tot_vel_indep_acc_;

    std::conditional_t<Interactions::ext_vel_indep, Coord, Empty> ext_vel_indep_acc_;

    std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> ext_vel_dep_acc_;
  };

/*---------------------------------------------------------------------------*\
    Class Accelerations Implementation
\*---------------------------------------------------------------------------*/
}
#endif //SPACEHUB_ACCELERATIONS_H
