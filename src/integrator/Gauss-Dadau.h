
#ifndef GAUSSDADAU_H
#define GAUSSDADAU_H

#include "../core-computation.hpp"
#include "../dev-tools.hpp"
#include "integrator.hpp"
#include <array>
#include <vector>

namespace space::integrator {

  class RadauConsts {
  public:
    [[nodiscard]] static constexpr double step_sequence(size_t i) { return sub_steps_[i]; }

    [[nodiscard]] static constexpr double G_tab(size_t n, size_t j) { return G_coef_[n * (n + 1) / 2 + j]; }

    [[nodiscard]] static constexpr double vel_B_tab(size_t stage, size_t i) { return vel_coef_[stage][i]; }

    [[nodiscard]] static constexpr double pos_B_tab(size_t stage, size_t i) { return pos_coef_[stage][i]; }

    template<typename Tab>
    void transfer_G_to_B(Tab const &G, Tab &B);

  private:
    static constexpr double sub_steps_[8] = {
            5.62625605369221465e-02, 1.80240691736892365e-01, 3.52624717113169637e-01, 5.47153626330555383e-01,
            7.34210177215410532e-01, 8.85320946839095768e-01, 9.77520613561287502e-01
    };

    static constexpr double G_coef_[28] = {
            17.7738089140780001,
            5.54813671853721663, 8.06593864838188711,
            2.83587607864443884, 3.37424997696263552, 5.80100155926406203,
            1.82764026751759771, 2.03711183535858464, 2.72544221180822598, 5.14062410581093269,
            1.36200781606246958, 1.47504021756041165, 1.80515358014025140, 2.62064492638703525, 5.34597689987110987,
            1.12953387533678987, 1.20618766605844559, 1.41827826373473910, 1.87724249618680994, 2.95711601729045588,
            6.61766201370242153,
            1.02299632982348675, 1.08547219393864239, 1.25426462228187776, 1.60026654949081622, 2.32359830021969444,
            4.10997577834455837, 10.8460261902368476
    };


    static constexpr double G_to_B_coef_[21] = {
            1.27179030902686775e-03, -1.43653023637089149e-03, 1.95656540994722109e-03,
            -3.57589772925161718e-03, 1.01408028300636298e-02, -5.62625605369221488e-02,
            -3.87603579159067693e-02, 4.21585277212687057e-02, -5.47553868890686898e-02,
            9.35376952594620670e-02, -2.36503252273814524e-01,
            3.60962243452845999e-01, -3.60099596502056807e-01, 4.15881200082306890e-01,
            -5.89127969386984196e-01,
            -1.46688420840042699e+00, 1.25015071184069093e+00, -1.13628159571753962e+00,
            2.90613625930842900e+00, -1.87049177293294999e+00,
            -2.75581271977204567e+00
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


  /** Gauss radau stepping. see details in
   https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
   */


  template<typename RadauTab>
  void transfer_G_to_B(const RadauTab &G, RadauTab &B) {
    size_t size = G.size();
    for (size_t i = 0; i < size; ++i) {
      B[i][0] = c[0] * G[i][6] + c[1] * G[i][5] + c[2] * G[i][4] + c[3] * G[i][3] + c[4] * G[i][2] + c[5] * G[i][1] +
                G[i][0];
      B[i][1] = c[6] * G[i][6] + c[7] * G[i][5] + c[8] * G[i][4] + c[9] * G[i][3] + c[10] * G[i][2] + G[i][1];
      B[i][2] = c[11] * G[i][6] + c[12] * G[i][5] + c[13] * G[i][4] + c[14] * G[i][3] + G[i][2];
      B[i][3] = c[15] * G[i][6] + c[16] * G[i][5] + c[17] * G[i][4] + G[i][3];
      B[i][4] = c[18] * G[i][6] + c[19] * G[i][5] + G[i][4];
      B[i][5] = c[20] * G[i][6] + G[i][5];
      B[i][6] = G[i][6];
    }
  }

  template<typename Coord>
  class GaussDadau : public Integrator<GaussDadau<Coord>> {
  public:
    using Base = Integrator<GaussDadau<Coord>>;
    using IterTable = std::array<Coord, 7>;
    using Scalar = typename Coord::Scalar;

    static constexpr size_t order{15};
    static constexpr size_t final_point{7};

    template<typename ParticleSys>
    void impl_integrate(ParticleSys const &particles, Scalar step_size);

    template<typename Scalar>
    void predict_new_B(Scalar step_ratio);

  private:
    template<typename ParticleSys>
    void integrate_to(ParticleSys const &particles, Scalar step_size, size_t stage);

    template<typename ParticleSys>
    void calc_B_table(ParticleSys &particles, Scalar step_size);

    void calc_G_table(Coord const &a0, Coord const a, size_t stage);

    void calc_vel_increament(Coord &dvel, Coord const &acc0, size_t stage);

    void calc_pos_increment(Coord &dpos, Coord const &vel0, Coord const &acc0, Scalar step_size, size_t stage);

  private:
    IterTable b_tab_;
    IterTable g_tab_;
    Coord acceleration0_;
    Coord acceleration_;
    Coord pos_increment_;
    Coord vel_increment_;
    std::vector<Scalar> input_;
  };

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::impl_integrate(const ParticleSys &particles, Scalar step_size) {
    particles.to_linear_container(input_);
    particles.evaluate_acc(acceleration0_);

    calc_B_table(particles, step_size);

    integrate_to(particles, step_size, final_point);
  }

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::integrate_to(const ParticleSys &particles, Scalar step_size, size_t stage) {
    calc_vel_increament(vel_increment_, acceleration0_, stage);
    calc_pos_increment(pos_increment_, particles.vel(), acceleration0_, step_size, stage);

    particles.advance_vel(step_size, vel_increment_);
    particles.advance_pos(step_size, pos_increment_);
    particles.advance_time(step_size);
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_vel_increament(Coord &dvel, const Coord &acc0, size_t stage) {
    calc::coord_scale(dvel, b_tab_[6], RadauConsts::vel_B_tab(stage, 7));
    calc::coord_advance(dvel, b_tab_[5], RadauConsts::vel_B_tab(stage, 6));
    calc::coord_advance(dvel, b_tab_[4], RadauConsts::vel_B_tab(stage, 5));
    calc::coord_advance(dvel, b_tab_[3], RadauConsts::vel_B_tab(stage, 4));
    calc::coord_advance(dvel, b_tab_[2], RadauConsts::vel_B_tab(stage, 3));
    calc::coord_advance(dvel, b_tab_[1], RadauConsts::vel_B_tab(stage, 2));
    calc::coord_advance(dvel, b_tab_[0], RadauConsts::vel_B_tab(stage, 1));
    calc::coord_advance(dvel, acc0, RadauConsts::vel_B_tab(stage, 0));
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_pos_increment(Coord &dpos, const Coord &vel0, const Coord &acc0, Scalar step_size,
                                             size_t stage) {
    calc::coord_scale(dpos, b_tab_[6], RadauConsts::pos_B_tab(stage, 8));
    calc::coord_advance(dpos, b_tab_[5], RadauConsts::pos_B_tab(stage, 7));
    calc::coord_advance(dpos, b_tab_[4], RadauConsts::pos_B_tab(stage, 6));
    calc::coord_advance(dpos, b_tab_[3], RadauConsts::pos_B_tab(stage, 5));
    calc::coord_advance(dpos, b_tab_[2], RadauConsts::pos_B_tab(stage, 4));
    calc::coord_advance(dpos, b_tab_[1], RadauConsts::pos_B_tab(stage, 3));
    calc::coord_advance(dpos, b_tab_[0], RadauConsts::pos_B_tab(stage, 2));
    calc::coord_advance(dpos, acc0, RadauConsts::pos_B_tab(stage, 1));
    calc::coord_scale(dpos, dpos, step_size);
    calc::coord_advance(dpos, vel0, RadauConsts::pos_B_tab(stage, 0));
  }

  template<typename Coord>
  void GaussDadau<Coord>::calc_G_table(Coord const &acc0, Coord const acc, size_t stage) {

    calc::coord_sub(g_tab_[stage], acc, acc0);

    calc::coord_scale(g_tab_[stage], g_tab_[stage], RadauConsts::G_tab(stage, 0));

    for (size_t j = 0; j < stage; ++j)
      calc::coord_advance(g_tab_[stage], g_tab_[j], -RadauConsts::G_tab(stage, j + 1));

  }

  template<typename Coord>
  template<typename ParticleSys>
  void GaussDadau<Coord>::calc_B_table(ParticleSys &particles, Scalar step_size) {
    integrate_to(particles, step_size, 0);
    particles.evaluate_acc(acceleration_);
    calc_G_table(acceleration0_, acceleration_, 0);

    for (size_t i = 1; i < final_point; ++i) {
      particles.load_from_linear_container(input_);
      integrate_to(particles, step_size, i);
      particles.evaluate_acc(acceleration_);
      calc_G_table(acceleration0_, acceleration_, i);
    }
    Radau::transfer_G_to_B(g_tab_, b_tab_);
    particles.load_from_linear_container(input_);
  }


/** @brief Gauss-Dadau integrator *//*
  template<typename ParticSys>
  class GaussDadau {
  public:
    *//* Typedef *//*
    SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
    using RadauArray = std::array<Vector, 7>;
    using RadauTab   = Container<RadauArray>;
    *//* Typedef *//*

    *//*Template parameter check*//*
    *//*Template parameter check*//*

    *//** @brief Order of the integrator*//*
    static const int order{15};
    static const size_t finalPoint{7};

    *//** @brief Interface to integrate particle system
     *
     *  This function integrate the particle system for one step with Gauss-Radau stepping.
     *  @param particles  Particle system need to be integrated.
     *  @param stepLength Step size for integration.
     *//*
    void integrate(ParticSys &particles, Scalar stepLength) {
      checkTabVolume(particles.particleNumber());
      calcuBTab(particles, stepLength);
      evaluateSystemAt(particles, stepLength, finalPoint);
    }

    *//**
     *
     * @param particles
     * @param stepLength
     *//*
    void calcuBTab(const ParticSys &particles, Scalar stepLength) {
      for (size_t i = 0; i < finalPoint; ++i) {
        localSystem_ = particles;
        evaluateSystemAt(localSystem_, stepLength, i);
        calcuGTab(particles.acc(), localSystem_.acc(), i);
      }
      Radau::transfer_G_to_B(Gtab_, Btab_);
    }

    *//**
     *
     * @param particleSys
     * @param stepLength
     * @param index
     *//*
    void evaluateSystemAt(ParticSys &particleSys, Scalar stepLength, size_t index) {
      VectorArray dpos;
      VectorArray dvel;
      evaluateVelIncrement(dvel, particleSys.acc(), index);
      evaluatePosIncrement(dpos, particleSys.vel(), particleSys.acc(), stepLength, index);
      particleSys.advanceVel(dvel, stepLength);
      particleSys.advancePos(dpos, stepLength);
      particleSys.advanceTime(stepLength);
      particleSys.evaluateAcc();
    }

    *//**
     *
     * @return
     *//*
    inline const RadauTab &getBTab() const {
      return Btab_;
    }

    *//**
     *
     * @return
     *//*
    inline const VectorArray &localAcc() const {
      return localSystem_.acc();
    }

    void checkTabVolume(size_t particleNum) {
      if constexpr (ParticSys::arraySize == space::DYNAMICAL) {

        if (particleNum != particleNumber_) {

          Btab_.resize(particleNum);
          Gtab_.resize(particleNum);
          oldBtab_.resize(particleNum);

          //once particles number changes, old b Value should be droped away.
          for (size_t i = 0; i < particleNum; ++i) {
            for (size_t j = 0; j < finalPoint; ++j) {
              Btab_[i][j].setZero();
              oldBtab_[i][j].setZero();
            }
          }
          particleNumber_ = particleNum;
        }
      }
    }

    *//**
     *
     * @param Q1
     *//*
    void predictNewB(Scalar Q1) {
      Scalar Q2 = Q1 * Q1;
      Scalar Q3 = Q2 * Q1;
      Scalar Q4 = Q2 * Q2;
      Scalar Q5 = Q3 * Q2;
      Scalar Q6 = Q3 * Q3;
      Scalar Q7 = Q4 * Q3;
      size_t size = Btab_.size();
      RadauArray dB;
      for (size_t i = 0; i < size; ++i) {

        for (size_t j = 0; j < finalPoint; ++j)
          dB[j] = Btab_[i][j] - oldBtab_[i][j];

        oldBtab_[i][0] = Q1 * (7 * Btab_[i][6] + 6 * Btab_[i][5] + 5 * Btab_[i][4] + 4 * Btab_[i][3] + 3 * Btab_[i][2] +
                               2 * Btab_[i][1] + Btab_[i][0]);
        oldBtab_[i][1] = Q2 *
                         (21 * Btab_[i][6] + 15 * Btab_[i][5] + 10 * Btab_[i][4] + 6 * Btab_[i][3] + 3 * Btab_[i][2] +
                          Btab_[i][1]);
        oldBtab_[i][2] = Q3 * (35 * Btab_[i][6] + 20 * Btab_[i][5] + 10 * Btab_[i][4] + 4 * Btab_[i][3] + Btab_[i][2]);
        oldBtab_[i][3] = Q4 * (35 * Btab_[i][6] + 15 * Btab_[i][5] + 5 * Btab_[i][4] + Btab_[i][3]);
        oldBtab_[i][4] = Q5 * (21 * Btab_[i][6] + 6 * Btab_[i][5] + Btab_[i][4]);
        oldBtab_[i][5] = Q6 * (7 * Btab_[i][6] + Btab_[i][5]);
        oldBtab_[i][6] = Q7 * Btab_[i][6];

        for (size_t j = 0; j < 7; ++j)
          Btab_[i][j] = oldBtab_[i][j] + dB[j];
      }
    }

  private:
    void evaluatePosIncrement(VectorArray &dpos, const VectorArray &vel, const VectorArray &acc, Scalar stepLength,
                              size_t iter) {
      using namespace Radau;
      size_t size = vel.size();
      dpos.resize(size);

      for (size_t i = 0; i < size; ++i) {
        dpos[i] = vel[i] * PC[iter][0] + (acc[i] * PC[iter][1] + Btab_[i][0] * PC[iter][2] + Btab_[i][1] * PC[iter][3] +
                                          Btab_[i][2] * PC[iter][4]
                                          + Btab_[i][3] * PC[iter][5] + Btab_[i][4] * PC[iter][6] +
                                          Btab_[i][5] * PC[iter][7] + Btab_[i][6] * PC[iter][8]) * stepLength;
      }
    }

    void evaluateVelIncrement(VectorArray &dvel, const VectorArray &acc, size_t iter) {
      using namespace Radau;
      size_t size = acc.size();
      dvel.resize(size);
      for (size_t i = 0; i < size; ++i) {
        dvel[i] =
                acc[i] * VC[iter][0] + Btab_[i][0] * VC[iter][1] + Btab_[i][1] * VC[iter][2] + Btab_[i][2] * VC[iter][3]
                + Btab_[i][3] * VC[iter][4] + Btab_[i][4] * VC[iter][5] + Btab_[i][5] * VC[iter][6] +
                Btab_[i][6] * VC[iter][7];
      }
    }

    *//**
     *
     * @param a0
     * @param a
     * @param iter
     *//*
    void calcuGTab(const VectorArray &a0, const VectorArray &a, size_t iter) {
      size_t size = a0.size();
      size_t offset = iter * (iter + 1) / 2 + 1;

      for (size_t i = 0; i < size; ++i) {
        Gtab_[i][iter] = (a[i] - a0[i]) * Radau::gg[offset - 1];

        for (size_t j = 0; j < iter; ++j)
          Gtab_[i][iter] -= Gtab_[i][j] * Radau::gg[j + offset];
      }
    }

  private:
    RadauTab oldBtab_;
    RadauTab Btab_;
    RadauTab Gtab_;
    ParticSys localSystem_;
    size_t particleNumber_{ParticSys::arraySize};
  };
*/
}
#endif
