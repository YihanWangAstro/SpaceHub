//
// Created by yihan on 5/16/19.
//

#ifndef SPACEHUB_SECULAR_H
#define SPACEHUB_SECULAR_H
/*---------------------------------------------------------------------------*\
        Class secular Declaration
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
        Class secular Implementation
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
        Help function
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
        Help macros
\*---------------------------------------------------------------------------*/

#include <array>
#include <iostream>
#include <tuple>
#include "/home/yihan/Repositories/odeint-v2/include/boost/numeric/odeint.hpp"
#include "vector3d.h"

using namespace boost::numeric::odeint;

namespace secular {

constexpr double pi = 3.14159265358979323;
constexpr double G = 4 * pi * pi;

double t_k_quad(double m_in, double m_out, double a_in, double a_out, double c_out) {
  double rc = a_out * c_out;
  return 1.0 / sqrt(G) * sqrt(m_in) / m_out * rc * rc * rc * pow(a_in, -1.5);
}

struct Secular {
 public:
  // Typemember
  using Container = std::array<double, 12>;

  Secular(double m1_, double m2_, double m3_)
      : m1{m1_}, m2{m2_}, m3{m3_}, mu1{m1 * m2 / (m1 + m2)}, mu2{(m1 + m2) * m3 / (m1 + m2 + m3)} {
    a_in_coef = 1 / sqrt(G * (m1 + m2)) / mu1;
    a_out_coef = 1 / sqrt(G * (m1 + m2 + m3)) / mu2;
  }

  void operator(Container const &x, Container &dxdt) {
    auto [L1, L2, e1, e2] = map_to_vectors(x);

    double L_in = norm(L1);

    double L_out = norm(L2);

    Vec3d n1 = L1 / L_in;

    Vec3d n2 = L2 / L_out;

    double de1e1 = norm2(e1);

    double de2e2 = norm2(e2);

    double dn1n2 = dot(n1, n2);

    double de1n2 = dot(e1, n2);

    Vec3d cn1n2 = cross(n1, n2);

    Vec3d cn1e1 = cross(n1, e1);

    Vec3d ce1n2 = cross(e1, n2);

    Vec3d ce2n1 = cross(e2, n1);

    Vec3d ce1e2 = cross(e1, e2);

    Vec3d cn2e2 = cross(n2, e2);

    double c_in = sqrt(1 - de1e1);

    double c_out = sqrt(1 - de2e2);

    double coef = 0.75 / t_k_quad(m1 + m2, m3, calc_a_in(L_in, c_in), calc_a_out(L_out, c_out), c_out);

    Vec3d dL1dt = coef * L_in * ((c_in * c_in * dn1n2) * cn1n2 - (5 * de1n2) * ce1n2);

    Vec3d dL2dt = -dL1dt;

    Vec3d de1dt = coef * c_in * (dn1n2 * ce1n2 + 2 * cn1e1 - 5 * de1n2 * cn1n2);

    Vec3d de2dt = coef * (L_in / L_out) / c_out *
                  ((c_in * c_in * dn1n2) * ce2n1 + (5 * de1n2) * ce1e2 -
                   (0.5 - 3 * de1e1 + 12.5 * de1n2 * de1n2 - 2.5 * c_in * c_in * dn1n2 * dn1n2) * cn2e2);

    map_to_array(dxdt, dj1, dj2, de1, de2);
  }

 private:
  auto map_to_vectors(state_type const &x) {
    return std::make_tuple(Vec3d{x[0], x[1], x[2]}, Vec3d{x[3], x[4], x[5]}, Vec3d{x[6], x[7], x[8]},
                           Vec3d{x[9], x[10], x[11]});
  }

  void map_to_array(state_type &x, Vec3d const &L1, Vec3d const &L2, Vec3d const &e1, Vec3d const &e2) {
    x[0] = L1.x, x[1] = L1.y, x[2] = L1.z;
    x[0] = L2.x, x[1] = L2.y, x[2] = L2.z;
    x[0] = e1.x, x[1] = e1.y, x[2] = e1.z;
    x[0] = e2.x, x[1] = e2.y, x[2] = e2.z;
  }

  double calc_a_in(double L_in, double c_in) { return a_in_coef * L_in / c_in; }

  double calc_a_out(double L_out, double c_out) { return a_out_coef * L_out / c_out; }

 private:
  double m1;
  double m2;
  double m3;
  double mu1;
  double mu2;
  double a_in_coef;
  double a_out_coef;
}

}  // namespace secular

void single_thread_run(state_type const &ini_conditions, double start_t, double end_t) {
  double ini_dt = (start_t - end_t) * 1e-5;
  auto observer = [&](state_type const &state, double t) {};
  integrate(secular, ini_conditions, start_t, end_t, ini_dt, ovserver);
}

int main(int argc, char **argv) {
  state_type ini_conditions = {10.0, 1.0, 1.0};
  double start_time = 0;
  double end_time = 25.0;
  single_thread_run(ini_conditions, start_time, end_time);
}

#endif  // SPACEHUB_SECULAR_H
