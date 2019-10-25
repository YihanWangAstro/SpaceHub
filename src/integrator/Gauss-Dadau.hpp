
#ifndef GAUSSDADAU_H
#define GAUSSDADAU_H

#include <array>
#include <vector>
#include "../core-computation.hpp"
#include "../dev-tools.hpp"
#include "integrator.hpp"

namespace space::integrator {

/*---------------------------------------------------------------------------*\
     Class RadauConsts Declaration
\*---------------------------------------------------------------------------*/
/**
 * Constant parameters used in Gauss Radau integration
 */
class RadauConsts {
 public:
  /**
   * Get the sub-step of the Gauss Radau integration.
   * @param[in] i The index of step sequence.
   * @return The sub-step.
   */
  [[nodiscard]] inline static constexpr double step_sequence(size_t i) { return sub_steps_[i]; }

  [[nodiscard]] inline static constexpr double G_tab(size_t n, size_t j) { return G_coef_[n * (n + 1) / 2 + j]; }

  [[nodiscard]] inline static constexpr double R_tab(size_t n, size_t j) { return R_coef_[n * (n + 1) / 2 + j]; }

  [[nodiscard]] inline static constexpr double B_tab(size_t n, size_t j) { return B_coef_[n * (n + 1) / 2 + j]; }

  [[nodiscard]] inline static constexpr double G2B_tab(size_t n, size_t j) { return G2B_coef_[n * (n + 1) / 2 + j]; }

  [[nodiscard]] inline static constexpr double B2G_tab(size_t n, size_t j) { return B2G_coef_[n * (n + 1) / 2 + j]; }

  [[nodiscard]] inline static constexpr double vel_B_tab(size_t stage, size_t i) { return vel_coef_[stage][i]; }

  [[nodiscard]] inline static constexpr double pos_B_tab(size_t stage, size_t i) { return pos_coef_[stage][i]; }

  template <typename Tab>
  static void transfer_G_to_B(Tab const &G, Tab &B);

  template <typename Tab>
  static void transfer_B_to_G(Tab const &B, Tab &G);

 private:
  static constexpr double sub_steps_[8] = {
      5.626256053692214647E-2, 1.802406917368923650E-1, 3.526247171131696374E-1, 5.471536263305553830E-1,
      7.342101772154105315E-1, 8.853209468390957681E-1, 9.775206135612875019E-1, 1};
  static constexpr double G_coef_[28] = {
      1.777380891407800084E+1, 4.475093038455599220E+1, 8.065938648381886689E+0, 5.550952167492268626E+1,
      1.957402937770697068E+1, 5.801001559264061482E+0, 5.216250225615303735E+1, 2.854090226792991111E+1,
      1.401047393301603805E+1, 5.140624105810934229E+0, 5.080809109074474632E+1, 3.730381756371242117E+1,
      2.529003421032797098E+1, 1.400990723922950856E+1, 5.345976899871107514E+0, 7.098538034164876460E+1,
      6.284484413580181971E+1, 5.210204506663947562E+1, 3.673612322693265978E+1, 1.956919433773405089E+1,
      6.617662013702424487E+0, 2.308581652314266533E+2, 2.256686153226572719E+2, 2.078990291808557011E+2,
      1.657537217326802767E+2, 1.035788205317551178E+2, 4.457690493316414861E+1, 1.084602619023684468E+1};
  static constexpr double R_coef_[28] = {
      1.777380891407800084E+1, 5.548136718537216506E+0, 8.065938648381886689E+0, 2.835876078644438678E+0,
      3.374249976962635260E+0, 5.801001559264061482E+0, 1.827640267517597830E+0, 2.037111835358584783E+0,
      2.725442211808226284E+0, 5.140624105810934229E+0, 1.362007816062469497E+0, 1.475040217560411548E+0,
      1.805153580140251260E+0, 2.620644926387035081E+0, 5.345976899871107514E+0, 1.129533875336789903E+0,
      1.206187666058445617E+0, 1.418278263734739154E+0, 1.877242496186810097E+0, 2.957116017290455748E+0,
      6.617662013702424487E+0, 1.022996329823486746E+0, 1.085472193938642384E+0, 1.254264622281877766E+0,
      1.600266549490816261E+0, 2.323598300219694223E+0, 4.109975778344559086E+0, 1.084602619023684468E+1};
  static constexpr double B_coef_[28] = {1, 2, 1,  3,  3,  1, 4, 6, 4,  1,  5,  10, 10, 5,
                                         1, 6, 15, 20, 15, 6, 1, 7, 21, 35, 35, 21, 7,  1};
  static constexpr double G2B_coef_[28] = {1,
                                           -5.62625605369221488e-02,
                                           1,
                                           1.01408028300636298e-02,
                                           -2.36503252273814524e-01,
                                           1,
                                           -3.57589772925161718e-03,
                                           9.35376952594620670e-02,
                                           -5.89127969386984196e-01,
                                           1,
                                           1.95656540994722109e-03,
                                           -5.47553868890686898e-02,
                                           4.15881200082306890e-01,
                                           -1.13628159571753962e+00,
                                           1,
                                           -1.43653023637089149e-03,
                                           4.21585277212687057e-02,
                                           -3.60099596502056807e-01,
                                           1.25015071184069093e+00,
                                           -1.87049177293294999e+00,
                                           1,
                                           1.27179030902686775e-03,
                                           -3.87603579159067693e-02,
                                           3.60962243452845999e-01,
                                           -1.46688420840042699e+00,
                                           2.90613625930842900e+00,
                                           -2.75581271977204567e+00,
                                           1};
  static constexpr double B2G_coef_[28] = {1,
                                           0.0562625605369221464656522,
                                           1,
                                           0.0031654757181708292499905,
                                           0.2365032522738145114532321,
                                           1,
                                           0.0001780977692217433881125,
                                           0.0457929855060279188954539,
                                           0.5891279693869841488271399,
                                           1,
                                           0.0000100202365223291272096,
                                           0.0084318571535257015445000,
                                           0.2535340690545692665214616,
                                           1.1362815957175395318285885,
                                           1,
                                           0.0000005637641639318207610,
                                           0.0015297840025004658189490,
                                           0.0978342365324440053653648,
                                           0.8752546646840910912297246,
                                           1.8704917729329500633517991,
                                           1,
                                           0.0000000317188154017613665,
                                           0.0002762930909826476593130,
                                           0.0360285539837364596003871,
                                           0.5767330002770787313544596,
                                           2.2485887607691597933926895,
                                           2.7558127197720458314421588,
                                           1};
  static constexpr double vel_coef_[8][8] = {
      {5.626256053692214647E-2, 1.582737859085414625E-3, 5.936592307391446270E-5, 2.505059130582281802E-6,
       1.127528327863641522E-7, 5.286469233626894413E-9, 2.549402531001510623E-10, 1.255064249542731518E-11},
      {1.802406917368923650E-1, 1.624335347889672983E-2, 1.951808844775469108E-3, 2.638465322403864957E-4,
       3.804470518671002818E-5, 5.714336649815626897E-6, 8.828194204973525680E-7, 1.392299851505546210E-7},
      {3.526247171131696374E-1, 6.217209555957145586E-2, 1.461561173935122318E-2, 3.865369466268484739E-3,
       1.090419851664646350E-3, 3.204241597731919312E-4, 9.684812459678303561E-5, 2.988216222152142116E-5},
      {5.471536263305553830E-1, 1.496885454033385145E-1, 5.460175362505510343E-2, 2.240666062496733591E-2,
       9.807908491927157006E-3, 4.472027248396827193E-3, 2.097330793722325170E-3, 1.004116880724923208E-3},
      {7.342101772154105315E-1, 2.695322921633422690E-1, 1.319289013296822230E-1, 7.264765651882529632E-2,
       4.267091901357679596E-2, 2.610785250908554176E-2, 1.643027230063671095E-2, 1.055536399535443914E-2},
      {8.853209468390957681E-1, 3.918965894560365175E-1, 2.313028397601537632E-1, 1.535829368272732326E-1,
       1.087761528402004554E-1, 8.025150552166705291E-2, 6.089857616031874067E-2, 4.717543696898040110E-2},
      {9.775206135612875019E-1, 4.777732749686179876E-1, 3.113554832603394526E-1, 2.282673022742386513E-1,
       1.785087947000749155E-1, 1.454133554344192817E-1, 1.218381877922220999E-1, 1.042119225751172763E-1},
      {1.000000000000000000E+0, 5.000000000000000000E-1, 3.333333333333333333E-1, 2.500000000000000000E-1,
       2.000000000000000000E-1, 1.666666666666666667E-1, 1.428571428571428571E-1, 1.250000000000000000E-1}};
  static constexpr double pos_coef_[8][9] = {
      {5.626256053692214647E-2, 1.582737859085414625E-3, 2.968296153695723135E-5, 8.350197101940939341E-7,
       2.818820819659103805E-8, 1.057293846725378883E-9, 4.249004218335851038E-11, 1.792948927918187883E-12,
       7.845903146402743936E-14},
      {1.802406917368923650E-1, 1.624335347889672983E-2, 9.759044223877345542E-4, 8.794884408012883189E-5,
       9.511176296677507045E-6, 1.142867329963125379E-6, 1.471365700828920947E-7, 1.988999787865066014E-8,
       2.788323203783690773E-9},
      {3.526247171131696374E-1, 6.217209555957145586E-2, 7.307805869675611591E-3, 1.288456488756161580E-3,
       2.726049629161615876E-4, 6.408483195463838623E-5, 1.614135409946383927E-5, 4.268880317360203023E-6,
       1.170798777788203990E-6},
      {5.471536263305553830E-1, 1.496885454033385145E-1, 2.730087681252755172E-2, 7.468886874989111969E-3,
       2.451977122981789251E-3, 8.944054496793654386E-4, 3.495551322870541951E-4, 1.434452686749890297E-4,
       6.104513250537416467E-5},
      {7.342101772154105315E-1, 2.695322921633422690E-1, 6.596445066484111152E-2, 2.421588550627509877E-2,
       1.066772975339419899E-2, 5.221570501817108351E-3, 2.738378716772785158E-3, 1.507909142193491306E-3,
       8.610950744002607232E-4},
      {8.853209468390957681E-1, 3.918965894560365175E-1, 1.156514198800768816E-1, 5.119431227575774419E-2,
       2.719403821005011384E-2, 1.605030110433341058E-2, 1.014976269338645678E-2, 6.739348138425771586E-3,
       4.640600280547312320E-3},
      {9.775206135612875019E-1, 4.777732749686179876E-1, 1.556777416301697263E-1, 7.608910075807955043E-2,
       4.462719867501872888E-2, 2.908267108688385633E-2, 2.030636463203701665E-2, 1.488741751073103947E-2,
       1.131881138844778091E-2},
      {1.000000000000000000E+0, 5.000000000000000000E-1, 1.666666666666666667E-1, 8.333333333333333333E-2,
       5.000000000000000000E-2, 3.333333333333333333E-2, 2.380952380952380952E-2, 1.785714285714285714E-2,
       1.388888888888888889E-2}};
};

/*---------------------------------------------------------------------------*\
     Class GaussRadau Declaration
\*---------------------------------------------------------------------------*/
/**
 * Gauss Radau stepping method. See details in
 * https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
 *
 * @tparam Coord
 */
template <typename Coord>
class GaussDadau : public Integrator<GaussDadau<Coord>> {
 public:
  using Base = Integrator<GaussDadau<Coord>>;
  using IterTable = std::array<Coord, 7>;
  using Scalar = typename Coord::Scalar;

  static constexpr size_t order{15};
  static constexpr size_t final_point{7};

  SPACEHUB_READ_ACCESSOR(auto, b_tab, b_tab_);

  SPACEHUB_READ_ACCESSOR(auto, init_acc, acceleration0_);

  SPACEHUB_READ_ACCESSOR(auto, last_acc, acceleration_);

  /**
   * Calculate the b table
   *
   * @tparam ParticleSys
   * @param particles
   * @param step_size
   */
  template <typename ParticleSys>
  void calc_B_table(ParticleSys &particles, Scalar step_size);

  void predict_new_B(Scalar step_ratio);

  template <typename ParticleSys>
  void integrate_to(ParticleSys &particles, Scalar step_size, size_t stage);

  void check_particle_size(size_t particle_num);

  CRTP_IMPL :

      template <typename ParticleSys>
      void
      impl_integrate(ParticleSys &particles, Scalar step_size);

 private:
  void calc_G_table(Coord const &acc0, Coord const &acc, size_t stage);

  void update_B_table(Coord const &acc0, Coord const &acc, size_t stage);

  void calc_vel_increment(Coord &dvel, Coord const &acc0, size_t stage);

  void calc_pos_increment(Coord &dpos, Coord const &vel0, Coord const &acc0, Scalar step_size, size_t stage);

 private:
  IterTable b_tab_;
  IterTable old_b_tab_;
  IterTable db_tab_;
  IterTable g_tab_;
  IterTable old_g_tab_;
  IterTable dg_tab_;

  Coord acceleration0_;
  Coord acceleration_;
  Coord pos_increment_;
  Coord vel_increment_;
  std::vector<Scalar> input_;
  size_t particle_num_{0};
};

/*---------------------------------------------------------------------------*\
     Class RadauConsts Implementation
\*---------------------------------------------------------------------------*/

template <typename Tab>
void RadauConsts::transfer_G_to_B(const Tab &G, Tab &B) {
  for (size_t stage = 0; stage < 7; ++stage) {
    calc::coord_scale(B[stage], G[6], G2B_tab(6, stage));
    for (size_t j = 6; j > stage; --j) {
      calc::coord_advance(B[stage], G[j - 1], G2B_tab(j - 1, stage));
    }
  }
}

template <typename Tab>
void RadauConsts::transfer_B_to_G(const Tab &B, Tab &G) {
  for (size_t stage = 0; stage < 7; ++stage) {
    calc::coord_scale(G[stage], B[6], B2G_tab(6, stage));
    for (size_t j = 6; j > stage; --j) {
      calc::coord_advance(G[stage], B[j - 1], B2G_tab(j - 1, stage));
    }
  }
}

/*---------------------------------------------------------------------------*\
     Class RadauConsts Implementation
\*---------------------------------------------------------------------------*/
template <typename Coord>
void GaussDadau<Coord>::check_particle_size(size_t particle_num) {
  if (particle_num_ != particle_num) {
    particle_num_ = particle_num;
    acceleration0_.resize(particle_num_);
    acceleration_.resize(particle_num_);
    pos_increment_.resize(particle_num_);
    vel_increment_.resize(particle_num_);
    for (auto &b : b_tab_) {
      b.resize(particle_num_);
      b.set_zero();
    }
    for (auto &old_b : old_b_tab_) {
      old_b.resize(particle_num_);
      old_b.set_zero();
    }
    for (auto &db : db_tab_) {
      db.resize(particle_num_);
      db.set_zero();
    }
    for (auto &g : g_tab_) {
      g.resize(particle_num_);
      g.set_zero();
    }
    for (auto &g : old_g_tab_) {
      g.resize(particle_num_);
      g.set_zero();
    }
    for (auto &g : dg_tab_) {
      g.resize(particle_num_);
      g.set_zero();
    }
  }
}

template <typename Coord>
template <typename ParticleSys>
void GaussDadau<Coord>::impl_integrate(ParticleSys &particles, Scalar step_size) {
  calc_B_table(particles, step_size);
  integrate_to(particles, step_size, final_point);
}

template <typename Coord>
template <typename ParticleSys>
void GaussDadau<Coord>::integrate_to(ParticleSys &particles, Scalar step_size, size_t stage) {
  calc_vel_increment(vel_increment_, acceleration0_, stage);
  calc_pos_increment(pos_increment_, particles.vel(), acceleration0_, step_size, stage);
  particles.advance_vel(step_size, vel_increment_);
  particles.advance_pos(step_size, pos_increment_);
  particles.advance_time(step_size);
}

template <typename Coord>
void GaussDadau<Coord>::calc_vel_increment(Coord &dvel, const Coord &acc0, size_t stage) {
  calc::coord_scale(dvel, b_tab_[6], RadauConsts::vel_B_tab(stage, 7));
  for (size_t i = 6; i > 0; --i) {
    calc::coord_advance(dvel, b_tab_[i - 1], RadauConsts::vel_B_tab(stage, i));
  }
  calc::coord_advance(dvel, acc0, RadauConsts::vel_B_tab(stage, 0));
}

template <typename Coord>
void GaussDadau<Coord>::calc_pos_increment(Coord &dpos, const Coord &vel0, const Coord &acc0, Scalar step_size,
                                           size_t stage) {
  calc::coord_scale(dpos, b_tab_[6], RadauConsts::pos_B_tab(stage, 8) * step_size);
  for (size_t i = 7; i > 1; --i) {
    calc::coord_advance(dpos, b_tab_[i - 2], RadauConsts::pos_B_tab(stage, i) * step_size);
  }
  calc::coord_advance(dpos, acc0, RadauConsts::pos_B_tab(stage, 1) * step_size);
  calc::coord_advance(dpos, vel0, RadauConsts::pos_B_tab(stage, 0));
}

template <typename Coord>
void GaussDadau<Coord>::calc_G_table(Coord const &acc0, Coord const &acc, size_t stage) {
  calc::coord_sub(g_tab_[stage], acc, acc0);
  calc::coord_scale(g_tab_[stage], g_tab_[stage], RadauConsts::G_tab(stage, 0));
  for (size_t j = 0; j < stage; ++j) {
    calc::coord_advance(g_tab_[stage], g_tab_[j], -RadauConsts::G_tab(stage, j + 1));
  }
}

template <typename Coord>
void GaussDadau<Coord>::update_B_table(const Coord &acc0, const Coord &acc, size_t stage) {
  calc::coord_sub(g_tab_[stage], acc, acc0);
  calc::coord_scale(g_tab_[stage], g_tab_[stage], RadauConsts::G_tab(stage, 0));
  for (size_t j = 0; j < stage; ++j) {
    calc::coord_advance(g_tab_[stage], old_g_tab_[j], -RadauConsts::G_tab(stage, j + 1));
  }

  calc::coord_sub(dg_tab_[stage], g_tab_[stage], old_g_tab_[stage]);

  for (size_t i = 0; i <= stage; ++i) calc::coord_advance(b_tab_[i], dg_tab_[stage], RadauConsts::G2B_tab(stage, i));

  swap(old_g_tab_[stage], g_tab_[stage]);
}

template <typename Coord>
template <typename ParticleSys>
void GaussDadau<Coord>::calc_B_table(ParticleSys &particles, Scalar step_size) {
  particles.to_linear_container(input_);
  check_particle_size(particles.number());
  particles.evaluate_acc(acceleration0_);
  integrate_to(particles, step_size, 0);
  particles.evaluate_acc(acceleration_);
  update_B_table(acceleration0_, acceleration_, 0);
  for (size_t i = 1; i < final_point; ++i) {
    particles.load_from_linear_container(input_);
    integrate_to(particles, step_size, i);
    particles.evaluate_acc(acceleration_);
    update_B_table(acceleration0_, acceleration_, i);
  }
  particles.load_from_linear_container(input_);
}

template <typename Coord>
void GaussDadau<Coord>::predict_new_B(Scalar step_ratio) {
  std::array<Scalar, final_point> Q;
  Q[0] = step_ratio;
  Q[1] = Q[0] * Q[0];
  Q[2] = Q[1] * Q[0];
  Q[3] = Q[1] * Q[1];
  Q[4] = Q[2] * Q[1];
  Q[5] = Q[2] * Q[2];
  Q[6] = Q[3] * Q[2];

  for (size_t j = 0; j < final_point; ++j) {
    calc::coord_sub(db_tab_[j], b_tab_[j], old_b_tab_[j]);
  }

  for (size_t stage = 0; stage < final_point; ++stage) {
    calc::coord_scale(old_b_tab_[stage], b_tab_[6], Q[stage] * RadauConsts::B_tab(6, stage));
    for (size_t j = 6; j > stage; --j) {
      calc::coord_advance(old_b_tab_[stage], b_tab_[j - 1], Q[stage] * RadauConsts::B_tab(j - 1, stage));
    }
  }

  for (size_t j = 0; j < final_point; ++j) {
    calc::coord_add(b_tab_[j], old_b_tab_[j], db_tab_[j]);
  }

  RadauConsts::transfer_B_to_G(b_tab_, old_g_tab_);
}
}  // namespace space::integrator
#endif
