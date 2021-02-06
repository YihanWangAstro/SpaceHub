
#pragma once

#include <array>
#include <vector>

#include "../core-computation.hpp"
#include "../dev-tools.hpp"

namespace space::integrator {

    /*---------------------------------------------------------------------------*\
         Class Radau Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Constant parameters used in Gauss Radau integration
     */
    class Radau {
       public:
        [[nodiscard]] inline static constexpr double h(size_t i) { return h_[i]; }

        [[nodiscard]] inline static constexpr double rs(size_t n, size_t j) { return rs_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double rr(size_t n, size_t j) { return rr_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double est_b(size_t n, size_t j) { return est_b_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double g2b(size_t n, size_t j) { return g2b_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double b2g(size_t n, size_t j) { return b2g_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double dy_tab(size_t stage, size_t i) { return dy_tab_[stage][i]; }

        template <typename Tab1, typename Tab2>
        static void transform_g2b(Tab1 const &G, Tab2 &B);

        template <typename Tab1, typename Tab2>
        static void transform_b2g(Tab1 const &B, Tab2 &G);

       private:
        static constexpr double h_[8] = {0.0562625605369221464656521910318, 0.180240691736892364987579942780,
                                         0.352624717113169637373907769648,  0.547153626330555383001448554766,
                                         0.734210177215410531523210605558,  0.885320946839095768090359771030,
                                         0.977520613561287501891174488626,  1.000000000000000000000000000000};
        static constexpr double est_b_[28] = {1.0, 2.0,  1.0,  3.0,  3.0,  1.0,  4.0,  6.0,  4.0,  1.0,
                                              5.0, 10.0, 10.0, 5.0,  1.0,  6.0,  15.0, 20.0, 15.0, 6.0,
                                              1.0, 7.0,  21.0, 35.0, 35.0, 21.0, 7.0,  1.0};
        static constexpr double g2b_[28] = {
            1.0000000000000000000000000000000,   -0.0562625605369221464656521910318, 1.0000000000000000000000000000000,
            0.0101408028300636299864818047860,   -0.236503252273814511453232133812,  1.0000000000000000000000000000000,
            -0.00357589772925161759493445889941, 0.0935376952594620658957484611455,  -0.589127969386984148827139903460,
            1.0000000000000000000000000000000,   0.00195656540994722107690056706032, -0.0547553868890686864408084294395,
            0.415881200082306861688621911192,    -1.13628159571753953182858845823,   1.0000000000000000000000000000000,
            -0.00143653023637089154244595529986, 0.0421585277212687077072973470356,  -0.360099596502056812289766461058,
            1.25015071184069102585054407511,     -1.87049177293295006335179906379,   1.0000000000000000000000000000000,
            0.00127179030902686774929431161484,  -0.0387603579159067703699046248206, 0.360962243452845983225339808035,
            -1.46688420840042696437015525831,    2.90613625930842930142379130730,    -2.75581271977204583144215883482,
            1.0000000000000000000000000000000};

        static constexpr double b2g_[28] = {1.0000000000000000000000000000000,    0.0562625605369221464656521910318,
                                            1.0000000000000000000000000000000,    0.00316547571817082924999048003940,
                                            0.236503252273814511453232133812,     1.0000000000000000000000000000000,
                                            0.000178097769221743388112527921974,  0.0457929855060279188954538730112,
                                            0.589127969386984148827139903460,     1.0000000000000000000000000000000,
                                            0.0000100202365223291272095672152244, 0.00843185715352570154449997416277,
                                            0.253534069054569266521461597106,     1.13628159571753953182858845823,
                                            1.0000000000000000000000000000000,    5.63764163931820761038385011543E-7,
                                            0.00152978400250046581894900795889,   0.0978342365324440053653648396422,
                                            0.875254664684091091229724588371,     1.87049177293295006335179906379,
                                            1.0000000000000000000000000000000,    3.17188154017613664758548178792E-8,
                                            0.000276293090982647659313022639369,  0.0360285539837364596003870741266,
                                            0.576733000277078731354459606135,     2.24858876076915979339268952601,
                                            2.75581271977204583144215883482,      1.0000000000000000000000000000000};

        static constexpr double rr_[28] = {
            0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.123978131199970218521927751748,
            0.352624717113169637373907769648,  0.296362156576247490908255578616, 0.172384025376277272386327826868,
            0.547153626330555383001448554766,  0.490891065793633236535796363735, 0.366912934593663018013868611986,
            0.194528909217385745627540785118,  0.734210177215410531523210605560, 0.677947616678488385057558414527,
            0.553969485478518166535630662777,  0.381585460102240894149302835910, 0.187056550884855148521762050792,
            0.885320946839095768090359771033,  0.829058386302173621624707579995, 0.705080255102203403102779828251,
            0.532696229725926130716452001382,  0.338167320508540385088911216264, 0.151110769623685236567149165472,
            0.977520613561287501891174488623,  0.921258053024365355425522297595, 0.797279921824395136903594545847,
            0.624895896448117864517266718977,  0.430366987230732118889725933859, 0.243310436345876970367963883068,
            0.0921996667221917338008147175960};

        static constexpr double rs_[28] = {
            17.7738089140780008407526623988, 44.7509303845559921986046250353, 8.06593864838188668853712230228,
            55.5095216749226862607771099218, 19.5740293777069706783363904372, 5.80100155926406148232868035040,
            52.1625022561530373477792209295, 28.5409022679299111073474153873, 14.0104739330160380493901680907,
            5.14062410581093422863632030242, 50.8080910907447463230817180713, 37.3038175637124211672269545084,
            25.2900342103279709751944749769, 14.0099072392295085641729477931, 5.34597689987110751412148951162,
            70.9853803416487645973553313623, 62.8448441358018197140240262907, 52.1020450666394756209272428984,
            36.7361232269326597782946191852, 19.5691943377340508924875261834, 6.61766201370242448744713000122,
            230.858165231426653338394745299, 225.668615322657271889813588120, 207.899029180855701108223308386,
            165.753721732680276715639824926, 103.578820531755117818385634006, 44.5769049331641486107191459286,
            10.8460261902368446847064289379};

        static constexpr double dy_tab_[8][8] = {
            {0.0562625605369221464656521910318, 0.00158273785908541462499524001970,
             0.0000593659230739144627041759739914, 0.00000250505913058228180239180380611,
             1.12752832786364152207677002309E-7, 5.28646923362689441264246964655E-9,
             2.54940253100151062270181599871E-10, 1.25506424954273151819358800414E-11},
            {0.180240691736892364987579942780, 0.0162433534788967298294907940929, 0.00195180884477546910833411890707,
             0.000263846532240386495674175398784, 0.0000380447051867100281790176392125,
             0.00000571433664981562689658801067631, 8.82819420497352567991757413913E-7,
             1.39229985150554621008628058320E-7},
            {0.352624717113169637373907769648, 0.0621720955595714558583705338679, 0.0146156117393512231817256964709,
             0.00386536946626848473895005296513, 0.00109041985166464635032388891393,
             0.000320424159773191931155202095438, 0.0000968481245967830356094210040723,
             0.0000298821622215214211602422306320},
            {0.547153626330555383001448554766, 0.149688545403338514457694770608, 0.0546017536250551034318076010137,
             0.0224066606249673359059495457779, 0.00980790849192715700573587464261, 0.00447202724839682719295985641144,
             0.00209733079372232517035152665694, 0.00100411688072492320816995569635},
            {0.734210177215410531523210605558, 0.269532292163342269000521386861, 0.131928901329682223034686180826,
             0.0726476565188252963150628540336, 0.0426709190135767959553564310090, 0.0261078525090855417556446504266,
             0.0164302723006367109483981032005, 0.0105553639953544391415383915456},
            {0.885320946839095768090359771030, 0.391896589456036517542394779982, 0.231302839760153763241855927822,
             0.153582936827273232576575659664, 0.108776152840200455366260812717, 0.0802515055216670529098759748592,
             0.0608985761603187406674123556033, 0.0471754369689804011020473412255},
            {0.977520613561287501891174488626, 0.477773274968617987575421375301, 0.311355483260339452636702056124,
             0.228267302274238651296621590670, 0.178508794700074915504569446223, 0.145413355434419281660632922534,
             0.121838187792222099875417308339, 0.104211922575117276292808124599},
            {1, 0.5, 0.333333333333333333333333333334, 0.250000000000000000000000000000,
             0.200000000000000000000000000000, 0.166666666666666666666666666667, 0.142857142857142857142857142857,
             0.125000000000000000000000000000}};
    };

    /*---------------------------------------------------------------------------*\
         Class GaussRadau Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Gauss Radau stepping method. See details in
     * https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
     *
     * @tparam TypeSystem
     */
    template <typename TypeSystem>
    class GaussRadau {
       public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        using State = AdScalarArray;
        using IterTable = std::array<ScalarArray, 7>;
        using IterStateTable = std::array<State, 7>;

        static constexpr size_t order{15};
        static constexpr size_t final_point{7};

        GaussRadau();

        SPACEHUB_READ_ACCESSOR(auto, b, b_);

        SPACEHUB_READ_ACCESSOR(auto, y_0, dydh0_);

        SPACEHUB_READ_ACCESSOR(auto, y_h, dydh_);

        SPACEHUB_READ_ACCESSOR(auto, diff_b6, dg_array_);  // after calc_b_table
        template <typename ParticleSys>
        void calc_b_table(ParticleSys &particles, Scalar step_size);

        void predict_new_b(Scalar step_ratio);

        template <typename ParticleSys>
        void integrate_to(ParticleSys &particles, Scalar step_size, size_t stage);

        template <typename ParticleSys>
        void integrate_at_end(ParticleSys &particles, Scalar step_size);

        void check_particle_size(size_t var_num);

        template <typename ParticleSys>
        void integrate(ParticleSys &particles, Scalar step_size);

       private:
        void update_b_table(ScalarArray const &y_init, ScalarArray const &y_now, size_t stage);

       private:
        IterTable b_;
        IterTable g_;
        IterTable old_b_;

        ScalarArray dydh0_{0};
        ScalarArray dydh_{0};
        ScalarArray tmp_array_{0};
        ScalarArray dg_array_{0};
        State tmp_state_{0};
        State input_{0};
        size_t var_num_{0};
    };

    /*---------------------------------------------------------------------------*\
         Class Radau Implementation
    \*---------------------------------------------------------------------------*/

    template <typename Tab1, typename Tab2>
    void Radau::transform_g2b(const Tab1 &G, Tab2 &B) {
        for (size_t stage = 0; stage < 7; ++stage) {
            calc::array_scale(B[stage], G[6], g2b(6, stage));
            for (size_t j = 6; j > stage; --j) {
                calc::array_advance(B[stage], G[j - 1], g2b(j - 1, stage));
            }
        }
    }

    template <typename Tab1, typename Tab2>
    void Radau::transform_b2g(const Tab1 &B, Tab2 &G) {
        for (size_t stage = 0; stage < 7; ++stage) {
            calc::array_scale(G[stage], B[6], b2g(6, stage));
            for (size_t j = 6; j > stage; --j) {
                calc::array_advance(G[stage], B[j - 1], b2g(j - 1, stage));
            }
        }
    }

    /*---------------------------------------------------------------------------*\
         Class Radau Implementation
    \*---------------------------------------------------------------------------*/

    template <typename TypeSystem>
    GaussRadau<TypeSystem>::GaussRadau() {
        auto init_iter_tab = [](auto &tab) {
            for (auto &t : tab) {
                t.resize(0);
            }
        };
        init_iter_tab(b_);
        init_iter_tab(old_b_);
        init_iter_tab(g_);
    }

    template <typename TypeSystem>
    void GaussRadau<TypeSystem>::check_particle_size(size_t var_num) {
        if (var_num_ != var_num) {
            var_num_ = var_num;

            resize_all(var_num_, dydh0_, dydh_, dg_array_, tmp_array_, tmp_state_);
            calc::set_arrays_zero(dydh0_, dydh_, dg_array_, tmp_array_, tmp_state_);

            auto set_iter_tab_0 = [](auto &tab, size_t num) {
                for (auto &t : tab) {
                    t.resize(num);
                    calc::array_set_zero(t);
                }
            };

            set_iter_tab_0(b_, var_num_);
            set_iter_tab_0(old_b_, var_num_);
            set_iter_tab_0(g_, var_num_);
        }
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussRadau<TypeSystem>::integrate(ParticleSys &particles, Scalar step_size) {
        calc_b_table(particles, step_size);
        integrate_at_end(particles, step_size);
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussRadau<TypeSystem>::integrate_to(ParticleSys &particles, Scalar step_size, size_t stage) {
        tmp_state_ = input_;
        auto h_n = Radau::h(stage);

        /*calc::array_scale(tmp_array_, b_[6], 7.0 * h_n / 8.0);
        for (size_t i = 6; i > 0; --i) {
            calc::array_scale_add(tmp_array_, tmp_array_, b_[i-1], static_cast<double>(i) * h_n / static_cast<double>(i
        + 1));
        }
        calc::array_advance(tmp_array_, dydh0_);
        calc::array_advance(tmp_state_, tmp_array_, step_size * h_n);*/
#pragma GCC ivdep
        for (size_t i = 0; i < var_num_; ++i) {
            tmp_state_[i] +=
                ((((((((b_[6][i] * (7.0 * h_n / 8.0)) + b_[5][i]) * (6.0 * h_n / 7.0) + b_[4][i]) * (5.0 * h_n / 6.0) +
                     b_[3][i]) *
                        (4.0 * h_n / 5.0) +
                    b_[2][i]) *
                       (3.0 * h_n / 4.0) +
                   b_[1][i]) *
                      (2.0 * h_n / 3.0) +
                  b_[0][i]) *
                     (1.0 * h_n / 2.0) +
                 dydh0_[i]) *
                (h_n * step_size);
        }
        particles.read_from_scalar_array(tmp_state_);
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussRadau<TypeSystem>::integrate_at_end(ParticleSys &particles, Scalar step_size) {
        tmp_state_ = input_;
        /* for (size_t i = 7; i > 0; --i) {
             calc::array_advance(tmp_state_, b_[i - 1], step_size / (i + 1));
         }
         calc::array_advance(tmp_state_, dydh0_, step_size);*/
#pragma GCC ivdep
        for (size_t i = 0; i < var_num_; ++i) {
            tmp_state_[i] += b_[6][i] * (step_size / 8.0);
            tmp_state_[i] += b_[5][i] * (step_size / 7.0);
            tmp_state_[i] += b_[4][i] * (step_size / 6.0);
            tmp_state_[i] += b_[3][i] * (step_size / 5.0);
            tmp_state_[i] += b_[2][i] * (step_size / 4.0);
            tmp_state_[i] += b_[1][i] * (step_size / 3.0);
            tmp_state_[i] += b_[0][i] * (step_size / 2.0);
            tmp_state_[i] += dydh0_[i] * step_size;
        }
        particles.read_from_scalar_array(tmp_state_);
    }

    template <typename TypeSystem>
    void GaussRadau<TypeSystem>::update_b_table(const ScalarArray &y_init, const ScalarArray &y_now, size_t stage) {
        /*calc::array_scale_sub(tmp_array_, y_now, y_init, Radau::rs(stage, 0));
        for (size_t j = 0; j < stage; ++j) {
            calc::array_advance(tmp_array_, g_[j], -Radau::rs(stage, j + 1));
        }
        calc::array_sub(dg_array_, tmp_array_, g_[stage]);

        std::swap(g_[stage], tmp_array_);

        for (size_t i = 0; i <= stage; ++i) {
            calc::array_advance(b_[i], dg_array_, Radau::g2b(stage, i));
        }*/
        switch (stage) {
            case 0:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[0][i];
                    g_[0][i] = (y_now[i] - y_init[i]) * Radau::rs(0, 0);
                    dg_array_[i] = g_[0][i] - tmp;
                }
                break;
            case 1:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[1][i];
                    g_[1][i] = (y_now[i] - y_init[i]) * Radau::rs(1, 0) - g_[0][i] * Radau::rs(1, 1);
                    dg_array_[i] = g_[1][i] - tmp;
                }
                break;
            case 2:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[2][i];
                    g_[2][i] = (y_now[i] - y_init[i]) * Radau::rs(2, 0) - g_[0][i] * Radau::rs(2, 1) -
                               g_[1][i] * Radau::rs(2, 2);
                    dg_array_[i] = g_[2][i] - tmp;
                }
                break;
            case 3:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[3][i];
                    g_[3][i] = (y_now[i] - y_init[i]) * Radau::rs(3, 0) - g_[0][i] * Radau::rs(3, 1) -
                               g_[1][i] * Radau::rs(3, 2) - g_[2][i] * Radau::rs(3, 3);
                    dg_array_[i] = g_[3][i] - tmp;
                }
                break;
            case 4:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[4][i];
                    g_[4][i] = (y_now[i] - y_init[i]) * Radau::rs(4, 0) - g_[0][i] * Radau::rs(4, 1) -
                               g_[1][i] * Radau::rs(4, 2) - g_[2][i] * Radau::rs(4, 3) - g_[3][i] * Radau::rs(4, 4);
                    dg_array_[i] = g_[4][i] - tmp;
                }
                break;
            case 5:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[5][i];
                    g_[5][i] = (y_now[i] - y_init[i]) * Radau::rs(5, 0) - g_[0][i] * Radau::rs(5, 1) -
                               g_[1][i] * Radau::rs(5, 2) - g_[2][i] * Radau::rs(5, 3) - g_[3][i] * Radau::rs(5, 4) -
                               g_[4][i] * Radau::rs(5, 5);
                    dg_array_[i] = g_[5][i] - tmp;
                }
                break;
            case 6:
#pragma GCC ivdep
                for (size_t i = 0; i < var_num_; ++i) {
                    auto tmp = g_[6][i];
                    g_[6][i] = (y_now[i] - y_init[i]) * Radau::rs(6, 0) - g_[0][i] * Radau::rs(6, 1) -
                               g_[1][i] * Radau::rs(6, 2) - g_[2][i] * Radau::rs(6, 3) - g_[3][i] * Radau::rs(6, 4) -
                               g_[4][i] * Radau::rs(6, 5) - g_[5][i] * Radau::rs(6, 6);
                    dg_array_[i] = g_[6][i] - tmp;
                }
                break;

            default:
                break;
        }
        for (size_t i = 0; i <= stage; ++i) {
            calc::array_advance(b_[i], dg_array_, Radau::g2b(stage, i));
        }
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussRadau<TypeSystem>::calc_b_table(ParticleSys &particles, Scalar step_size) {
        particles.write_to_scalar_array(input_);
        check_particle_size(input_.size());
        particles.evaluate_general_derivative(dydh0_);

        for (size_t i = 0; i < final_point; ++i) {
            integrate_to(particles, step_size, i);
            particles.evaluate_general_derivative(dydh_);
            update_b_table(dydh0_, dydh_, i);
            particles.read_from_scalar_array(input_);
        }
    }

    template <typename TypeSystem>
    void GaussRadau<TypeSystem>::predict_new_b(Scalar step_ratio) {
        std::array<Scalar, final_point> Q;
        Q[0] = step_ratio;
        Q[1] = Q[0] * Q[0];
        Q[2] = Q[1] * Q[0];
        Q[3] = Q[1] * Q[1];
        Q[4] = Q[2] * Q[1];
        Q[5] = Q[2] * Q[2];
        Q[6] = Q[3] * Q[2];

        for (size_t i = 0; i < final_point; ++i) {
            calc::array_sub(tmp_array_, b_[i], old_b_[i]);
            calc::array_scale(old_b_[i], b_[6], Radau::est_b(6, i));
            for (size_t j = 6; j > i; --j) {
                calc::array_advance(old_b_[i], b_[j - 1], Radau::est_b(j - 1, i));
            }
            calc::array_scale(old_b_[i], old_b_[i], Q[i]);
            calc::array_add(b_[i], old_b_[i], tmp_array_);
        }
        Radau::transform_b2g(b_, g_);
    }
}  // namespace space::integrator
