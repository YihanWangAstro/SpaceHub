
#ifndef GAUSSDADAU_H
#define GAUSSDADAU_H

#include "../core-computation.hpp"
#include "../dev-tools.hpp"
#include "integrator.hpp"
#include <array>
#include <vector>
/** Gauss radau stepping. see details in
   https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
   */
namespace space::integrator {
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

    [[nodiscard]] inline static constexpr double B_tab(size_t n, size_t j) { return B_coef_[n * (n + 1) / 2 + j]; }

    [[nodiscard]] inline static constexpr double G2B_tab(size_t n, size_t j) { return G2B_coef_[n * (n + 1) / 2 + j]; }

    [[nodiscard]] inline static constexpr double vel_B_tab(size_t stage, size_t i) { return vel_coef_[stage][i]; }

    [[nodiscard]] inline static constexpr double pos_B_tab(size_t stage, size_t i) { return pos_coef_[stage][i]; }

    template<typename Tab>
    static void transfer_G_to_B(Tab const &G, Tab &B) {
      for (size_t stage = 0; stage < 7; ++stage) {
        calc::coord_scale(B[stage], G[6], G2B_tab(6, stage));
        for (size_t j = 6; j > stage; --j) {
          calc::coord_advance(B[stage], G[j - 1], G2B_tab(j - 1, stage));
        }
      }
    }

  private:
    static constexpr double sub_steps_[8] = {
            5.62625605369221465e-02, 1.80240691736892365e-01, 3.52624717113169637e-01, 5.47153626330555383e-01,
            7.34210177215410532e-01, 8.85320946839095768e-01, 9.77520613561287502e-01
    };
    static constexpr double G_coef_[28] = {
            1.777380891407800010E+1,
            4.475093038455599560E+1, 8.065938648381887110E+0,
            5.550952167492269900E+1, 1.957402937770697400E+1, 5.801001559264062030E+0,
            5.216250225615300860E+1, 2.854090226792989740E+1, 1.401047393301603230E+1, 5.140624105810932690E+0,
            5.080809109074478230E+1, 3.730381756371244540E+1, 2.529003421032798570E+1, 1.400990723922951560E+1,
            5.345976899871109870E+0,
            7.098538034164872380E+1, 6.284484413580178570E+1, 5.210204506663944840E+1, 3.673612322693264200E+1,
            1.956919433773404300E+1, 6.617662013702421530E+0,
            2.308581652314266940E+2, 2.256686153226573090E+2, 2.078990291808557330E+2, 1.657537217326803030E+2,
            1.035788205317551370E+2, 4.457690493316415290E+1, 1.084602619023684760E+1
    };
    static constexpr double B_coef_[28] = {
            1,
            2, 1,
            3, 3, 1,
            4, 6, 4, 1,
            5, 10, 10, 5, 1,
            6, 15, 20, 15, 6, 1,
            7, 21, 35, 35, 21, 7, 1
    };
    static constexpr double G2B_coef_[28] = {
            1,
            -5.62625605369221488e-02, 1,
            1.01408028300636298e-02, -2.36503252273814524e-01, 1,
            -3.57589772925161718e-03, 9.35376952594620670e-02, -5.89127969386984196e-01, 1,
            1.95656540994722109e-03, -5.47553868890686898e-02, 4.15881200082306890e-01, -1.13628159571753962e+00, 1,
            -1.43653023637089149e-03, 4.21585277212687057e-02, -3.60099596502056807e-01, 1.25015071184069093e+00,
            -1.87049177293294999e+00, 1,
            1.27179030902686775e-03, -3.87603579159067693e-02, 3.60962243452845999e-01, -1.46688420840042699e+00,
            2.90613625930842900e+00, -2.75581271977204567e+00, 1
    };
    static constexpr double vel_coef_[8][8] = {
            {5.626256053692214880E-2, 1.582737859085414760E-3, 5.936592307391447020E-5, 2.505059130582282220E-6, 1.127528327863641760E-7, 5.286469233626895740E-9, 2.549402531001511370E-10, 1.255064249542731940E-11},
            {1.802406917368923610E-1, 1.624335347889672910E-2, 1.951808844775468980E-3, 2.638465322403864740E-4, 3.804470518671002400E-5, 5.714336649815626160E-6, 8.828194204973524360E-7,  1.392299851505545970E-7},
            {3.526247171131696170E-1, 6.217209555957144850E-2, 1.461561173935122060E-2, 3.865369466268483810E-3, 1.090419851664646030E-3, 3.204241597731918180E-4, 9.684812459678299580E-5,  2.988216222152140710E-5},
            {5.471536263305554200E-1, 1.496885454033385350E-1, 5.460175362505511470E-2, 2.240666062496734200E-2, 9.807908491927160340E-3, 4.472027248396829020E-3, 2.097330793722326170E-3,  1.004116880724923760E-3},
            {7.342101772154104870E-1, 2.695322921633422370E-1, 1.319289013296821990E-1, 7.264765651882527860E-2, 4.267091901357678300E-2, 2.610785250908553220E-2, 1.643027230063670390E-2,  1.055536399535443400E-2},
            {8.853209468390957900E-1, 3.918965894560365370E-1, 2.313028397601537810E-1, 1.535829368272732480E-1, 1.087761528402004690E-1, 8.025150552166706500E-2, 6.089857616031875140E-2,  4.717543696898041060E-2},
            {9.775206135612874990E-1, 4.777732749686179850E-1, 3.113554832603394500E-1, 2.282673022742386490E-1, 1.785087947000749130E-1, 1.454133554344192790E-1, 1.218381877922220970E-1,  1.042119225751172740E-1},
            {1.000000000000000000E+0, 5.000000000000000000E-1, 3.333333333333333340E-1, 2.500000000000000000E-1, 2.000000000000000000E-1, 1.666666666666666670E-1, 1.428571428571428570E-1,  1.250000000000000000E-1}
    };
    static constexpr double pos_coef_[8][9] = {
            {5.626256053692214880E-2, 1.582737859085414760E-3, 2.968296153695723500E-5, 8.350197101940940710E-7, 2.818820819659104380E-8, 1.057293846725379140E-9, 4.249004218335852240E-11, 1.792948927918188470E-12, 7.845903146402746830E-14},
            {1.802406917368923610E-1, 1.624335347889672910E-2, 9.759044223877344890E-4, 8.794884408012882400E-5, 9.511176296677506010E-6, 1.142867329963125230E-6, 1.471365700828920740E-7,  1.988999787865065700E-8,  2.788323203783690250E-9},
            {3.526247171131696170E-1, 6.217209555957144850E-2, 7.307805869675610290E-3, 1.288456488756161280E-3, 2.726049629161615060E-4, 6.408483195463836340E-5, 1.614135409946383260E-5,  4.268880317360200990E-6,  1.170798777788203370E-6},
            {5.471536263305554200E-1, 1.496885454033385350E-1, 2.730087681252755730E-2, 7.468886874989113980E-3, 2.451977122981790090E-3, 8.944054496793658040E-4, 3.495551322870543610E-4,  1.434452686749891080E-4,  6.104513250537420200E-5},
            {7.342101772154104870E-1, 2.695322921633422370E-1, 6.596445066484109930E-2, 2.421588550627509280E-2, 1.066772975339419570E-2, 5.221570501817106420E-3, 2.738378716772783980E-3,  1.507909142193490570E-3,  8.610950744002602490E-4},
            {8.853209468390957900E-1, 3.918965894560365370E-1, 1.156514198800768900E-1, 5.119431227575774920E-2, 2.719403821005011720E-2, 1.605030110433341290E-2, 1.014976269338645850E-2,  6.739348138425772970E-3,  4.640600280547313390E-3},
            {9.775206135612874990E-1, 4.777732749686179850E-1, 1.556777416301697250E-1, 7.608910075807954960E-2, 4.462719867501872830E-2, 2.908267108688385580E-2, 2.030636463203701620E-2,  1.488741751073103920E-2,  1.131881138844778050E-2},
            {1.000000000000000000E+0, 5.000000000000000000E-1, 1.666666666666666660E-1, 8.333333333333333300E-2, 4.999999999999999980E-2, 3.333333333333333320E-2, 2.380952380952380940E-2,  1.785714285714285700E-2,  1.388888888888888880E-2}
    };
  };

  template<typename Coord>
  class GaussDadau : public Integrator<GaussDadau<Coord>> {
  public:
    using Base = Integrator<GaussDadau<Coord>>;
    using IterTable = std::array<Coord, 7>;
    using Scalar = typename Coord::Scalar;

    static constexpr size_t order{15};
    static constexpr size_t final_point{7};

    SPACEHUB_READ_ACCESSOR(auto, b_tab, b_tab_);

    SPACEHUB_READ_ACCESSOR(auto, init_acc, acceleration0_);

    template<typename ParticleSys>
    void impl_integrate(ParticleSys &particles, Scalar step_size);

    template<typename ParticleSys>
    void calc_B_table(ParticleSys &particles, Scalar step_size);

    void predict_new_B(Scalar step_ratio);

    template<typename ParticleSys>
    void integrate_to(ParticleSys &particles, Scalar step_size, size_t stage);

    void check_particle_size(size_t particle_num);

  private:
    void calc_G_table(Coord const &acc0, Coord const &acc, size_t stage);

    void calc_vel_increment(Coord &dvel, Coord const &acc0, size_t stage);

    void calc_pos_increment(Coord &dpos, Coord const &vel0, Coord const &acc0, Scalar step_size, size_t stage);

  private:
    IterTable b_tab_;
    IterTable old_b_tab_;
    IterTable db_tab_;
    IterTable g_tab_;
    Coord acceleration0_;
    Coord acceleration_;
    Coord pos_increment_;
    Coord vel_increment_;
    std::vector<Scalar> input_;
    size_t particle_num_{0};
  };

  template<typename Coord>
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
    }
  }

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::impl_integrate(ParticleSys &particles, Scalar step_size) {
    calc_B_table(particles, step_size);
    integrate_to(particles, step_size, final_point);
  }

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::integrate_to(ParticleSys &particles, Scalar step_size, size_t stage) {
    calc_vel_increment(vel_increment_, acceleration0_, stage);
    calc_pos_increment(pos_increment_, particles.vel(), acceleration0_, step_size, stage);
    particles.advance_vel(step_size, vel_increment_);
    particles.advance_pos(step_size, pos_increment_);
    particles.advance_time(step_size);
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_vel_increment(Coord &dvel, const Coord &acc0, size_t stage) {
    calc::coord_scale(dvel, b_tab_[6], RadauConsts::vel_B_tab(stage, 7));
    for (size_t i = 6; i > 0; --i) {
      calc::coord_advance(dvel, b_tab_[i - 1], RadauConsts::vel_B_tab(stage, i));
    }
    calc::coord_advance(dvel, acc0, RadauConsts::vel_B_tab(stage, 0));
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_pos_increment(Coord &dpos, const Coord &vel0, const Coord &acc0, Scalar step_size,
                                             size_t stage) {
    calc::coord_scale(dpos, b_tab_[6], RadauConsts::pos_B_tab(stage, 8));
    for (size_t i = 7; i > 1; --i) {
      calc::coord_advance(dpos, b_tab_[i - 2], RadauConsts::pos_B_tab(stage, i));
    }
    calc::coord_advance(dpos, acc0, RadauConsts::pos_B_tab(stage, 1));
    calc::coord_scale(dpos, dpos, step_size);
    calc::coord_advance(dpos, vel0, RadauConsts::pos_B_tab(stage, 0));
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_G_table(Coord const &acc0, Coord const &acc, size_t stage) {
    calc::coord_sub(g_tab_[stage], acc, acc0);
    calc::coord_scale(g_tab_[stage], g_tab_[stage], RadauConsts::G_tab(stage, 0));
    for (size_t j = 0; j < stage; ++j) {
      calc::coord_advance(g_tab_[stage], g_tab_[j], -RadauConsts::G_tab(stage, j + 1));
    }
  }

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::calc_B_table(ParticleSys &particles, Scalar step_size) {
    particles.to_linear_container(input_);
    check_particle_size(particles.number());
    particles.evaluate_acc(acceleration0_);
    integrate_to(particles, step_size, 0);
    particles.evaluate_acc(acceleration_);
    calc_G_table(acceleration0_, acceleration_, 0);
    for (size_t i = 1; i < final_point; ++i) {
      particles.load_from_linear_container(input_);
      integrate_to(particles, step_size, i);
      particles.evaluate_acc(acceleration_);
      calc_G_table(acceleration0_, acceleration_, i);
    }
    RadauConsts::transfer_G_to_B(g_tab_, b_tab_);
    particles.load_from_linear_container(input_);
  }

  template<typename Coord>
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
  }
}
#endif
