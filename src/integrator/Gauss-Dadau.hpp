
#pragma once

#include <array>
#include <vector>

#include "../core-computation.hpp"
#include "../dev-tools.hpp"

namespace space::integrator {

    /*---------------------------------------------------------------------------*\
         Class RadauConsts Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * Constant parameters used in Gauss Radau integration
     */
    class RadauConsts {
       public:
        [[nodiscard]] inline static constexpr double step_sequence(size_t i) { return sub_steps_[i]; }

        [[nodiscard]] inline static constexpr double G_tab(size_t n, size_t j) { return G_coef_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double RR_tab(size_t n, size_t j) { return rr_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double B_tab(size_t n, size_t j) { return B_coef_[n * (n + 1) / 2 + j]; }

        [[nodiscard]] inline static constexpr double G2B_tab(size_t n, size_t j) {
            return G2B_coef_[n * (n + 1) / 2 + j];
        }

        [[nodiscard]] inline static constexpr double B2G_tab(size_t n, size_t j) {
            return B2G_coef_[n * (n + 1) / 2 + j];
        }

        [[nodiscard]] inline static constexpr double dy_tab(size_t stage, size_t i) { return dy_tab_[stage][i]; }

        template <typename Tab1, typename Tab2>
        static void transfer_G_to_B(Tab1 const &G, Tab2 &B);

        template <typename Tab1, typename Tab2>
        static void transfer_B_to_G(Tab1 const &B, Tab2 &G);

       private:
        static constexpr double sub_steps_[8] = {0.0562625605369221464656521910318, 0.180240691736892364987579942780,
                                                 0.352624717113169637373907769648, 0.547153626330555383001448554766,
                                                 0.734210177215410531523210605558, 0.885320946839095768090359771030,
                                                 0.977520613561287501891174488626, 1.000000000000000000000000000000};
        static constexpr double B_coef_[28] = {1, 2, 1, 3, 3, 1, 4, 6, 4, 1, 5, 10, 10, 5,
                                               1, 6, 15, 20, 15, 6, 1, 7, 21, 35, 35, 21, 7, 1};
        static constexpr double G2B_coef_[28] = {
                1.0000000000000000000000000, -0.0562625605369221464656522, 1.0000000000000000000000000,
                0.0101408028300636299864818, -0.2365032522738145114532321, 1.0000000000000000000000000,
                -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399,
                1.0000000000000000000000000, 0.0019565654099472210769006, -0.0547553868890686864408084,
                0.4158812000823068616886219, -1.1362815957175395318285885, 1.0000000000000000000000000,
                -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665,
                1.2501507118406910258505441, -1.8704917729329500633517991, 1.0000000000000000000000000,
                0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398,
                -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588,
                1.0000000000000000000000000};

        static constexpr double B2G_coef_[28] = {
                1.0000000000000000000000000, 0.0562625605369221464656522, 1.0000000000000000000000000,
                0.0031654757181708292499905, 0.2365032522738145114532321, 1.0000000000000000000000000,
                0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399,
                1.0000000000000000000000000, 0.0000100202365223291272096, 0.0084318571535257015445000,
                0.2535340690545692665214616, 1.1362815957175395318285885, 1.0000000000000000000000000,
                0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648,
                0.8752546646840910912297246, 1.8704917729329500633517991, 1.0000000000000000000000000,
                0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871,
                0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588,
                1.0000000000000000000000000};

        static constexpr double rr_[28] = {
                0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278,
                0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278,
                0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686,
                0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584,
                0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621,
                0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798,
                0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492,
                0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945,
                0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639,
                0.0921996667221917338008147};

        static constexpr double G_coef_[28] = {
                1.777380891407800084E+1, 4.475093038455599220E+1, 8.065938648381886689E+0, 5.550952167492268626E+1,
                1.957402937770697068E+1, 5.801001559264061482E+0, 5.216250225615303735E+1, 2.854090226792991111E+1,
                1.401047393301603805E+1, 5.140624105810934229E+0, 5.080809109074474632E+1, 3.730381756371242117E+1,
                2.529003421032797098E+1, 1.400990723922950856E+1, 5.345976899871107514E+0, 7.098538034164876460E+1,
                6.284484413580181971E+1, 5.210204506663947562E+1, 3.673612322693265978E+1, 1.956919433773405089E+1,
                6.617662013702424487E+0, 2.308581652314266533E+2, 2.256686153226572719E+2, 2.078990291808557011E+2,
                1.657537217326802767E+2, 1.035788205317551178E+2, 4.457690493316414861E+1, 1.084602619023684468E+1};

        static constexpr double dy_tab_[8][8] = {
                {5.626256053692214647E-2, 1.582737859085414625E-3, 5.936592307391446270E-5, 2.505059130582281802E-6,
                        1.127528327863641522E-7, 5.286469233626894413E-9, 2.549402531001510623E-10, 1.255064249542731518E-11},
                {1.802406917368923650E-1, 1.624335347889672983E-2, 1.951808844775469108E-3, 2.638465322403864957E-4,
                        3.804470518671002818E-5, 5.714336649815626897E-6, 8.828194204973525680E-7,  1.392299851505546210E-7},
                {3.526247171131696374E-1, 6.217209555957145586E-2, 1.461561173935122318E-2, 3.865369466268484739E-3,
                        1.090419851664646350E-3, 3.204241597731919312E-4, 9.684812459678303561E-5,  2.988216222152142116E-5},
                {5.471536263305553830E-1, 1.496885454033385145E-1, 5.460175362505510343E-2, 2.240666062496733591E-2,
                        9.807908491927157006E-3, 4.472027248396827193E-3, 2.097330793722325170E-3,  1.004116880724923208E-3},
                {7.342101772154105315E-1, 2.695322921633422690E-1, 1.319289013296822230E-1, 7.264765651882529632E-2,
                        4.267091901357679596E-2, 2.610785250908554176E-2, 1.643027230063671095E-2,  1.055536399535443914E-2},
                {8.853209468390957681E-1, 3.918965894560365175E-1, 2.313028397601537632E-1, 1.535829368272732326E-1,
                        1.087761528402004554E-1, 8.025150552166705291E-2, 6.089857616031874067E-2,  4.717543696898040110E-2},
                {9.775206135612875019E-1, 4.777732749686179876E-1, 3.113554832603394526E-1, 2.282673022742386513E-1,
                        1.785087947000749155E-1, 1.454133554344192817E-1, 1.218381877922220999E-1,  1.042119225751172763E-1},
                {1.000000000000000000E+0, 5.000000000000000000E-1, 3.333333333333333333E-1, 2.500000000000000000E-1,
                        2.000000000000000000E-1, 1.666666666666666667E-1, 1.428571428571428571E-1,  1.250000000000000000E-1}};
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
    class GaussDadau {
       public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        using State = AdScalarArray;
        using IterTable = std::array<ScalarArray, 7>;
        using IterStateTable = std::array<State, 7>;

        static constexpr size_t order{15};
        static constexpr size_t final_point{7};

        GaussDadau();

        SPACEHUB_READ_ACCESSOR(auto, b_tab, b_tab_);

        SPACEHUB_READ_ACCESSOR(auto, init_acc, acceleration0_);

        SPACEHUB_READ_ACCESSOR(auto, last_acc, acceleration_);

        SPACEHUB_READ_ACCESSOR(auto, diff_b6, dg_array_);  // calc_B_table
        template <typename ParticleSys>
        void calc_B_table(ParticleSys &particles, Scalar step_size);

        void predict_new_B(Scalar step_ratio);

        template <typename ParticleSys>
        void integrate_to(ParticleSys &particles, Scalar step_size, size_t stage);

        template <typename ParticleSys>
        void integrate_at_end(ParticleSys &particles, Scalar step_size);

        void check_particle_size(size_t var_num);

        template <typename ParticleSys>
        void integrate(ParticleSys &particles, Scalar step_size);

       private:
        void update_B_table(ScalarArray const &acc0, ScalarArray const &acc, size_t stage);

       private:
        IterTable b_tab_;
        IterTable old_b_tab_;
        IterTable g_tab_;

        ScalarArray acceleration0_{0};
        ScalarArray acceleration_{0};
        ScalarArray tmp_array_{0};
        ScalarArray dg_array_{0};
        State tmp_state_{0};
        State input_{0};
        size_t var_num_{0};
    };

    /*---------------------------------------------------------------------------*\
         Class RadauConsts Implementation
    \*---------------------------------------------------------------------------*/

    template <typename Tab1, typename Tab2>
    void RadauConsts::transfer_G_to_B(const Tab1 &G, Tab2 &B) {
        for (size_t stage = 0; stage < 7; ++stage) {
            calc::array_scale(B[stage], G[6], G2B_tab(6, stage));
            for (size_t j = 6; j > stage; --j) {
                calc::array_advance(B[stage], G[j - 1], G2B_tab(j - 1, stage));
            }
        }
    }

    template <typename Tab1, typename Tab2>
    void RadauConsts::transfer_B_to_G(const Tab1 &B, Tab2 &G) {
        for (size_t stage = 0; stage < 7; ++stage) {
            calc::array_scale(G[stage], B[6], B2G_tab(6, stage));
            for (size_t j = 6; j > stage; --j) {
                calc::array_advance(G[stage], B[j - 1], B2G_tab(j - 1, stage));
            }
        }
    }

    /*---------------------------------------------------------------------------*\
         Class RadauConsts Implementation
    \*---------------------------------------------------------------------------*/

    template <typename TypeSystem>
    GaussDadau<TypeSystem>::GaussDadau() {
        auto init_iter_tab = [](auto &tab) {
            for (auto &t : tab) {
                t.resize(0);
            }
        };
        init_iter_tab(b_tab_);
        init_iter_tab(old_b_tab_);
        init_iter_tab(g_tab_);
    }

    template <typename TypeSystem>
    void GaussDadau<TypeSystem>::check_particle_size(size_t var_num) {
        if (var_num_ != var_num) {
            var_num_ = var_num;

            resize_all(var_num_, acceleration0_, acceleration_, dg_array_, tmp_array_, tmp_state_);
            calc::set_arrays_zero(acceleration0_, acceleration_, dg_array_, tmp_array_, tmp_state_);

            auto set_iter_tab_0 = [](auto &tab, size_t num) {
                for (auto &t : tab) {
                    t.resize(num);
                    calc::array_set_zero(t);
                }
            };

            set_iter_tab_0(b_tab_, var_num_);
            set_iter_tab_0(old_b_tab_, var_num_);
            set_iter_tab_0(g_tab_, var_num_);
        }
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussDadau<TypeSystem>::integrate(ParticleSys &particles, Scalar step_size) {
        calc_B_table(particles, step_size);
        integrate_at_end(particles, step_size);
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussDadau<TypeSystem>::integrate_to(ParticleSys &particles, Scalar step_size, size_t stage) {
        tmp_state_ = input_;

        auto h_n = RadauConsts::step_sequence(stage);
        calc::array_scale(tmp_array_, b_tab_[6], 7 * h_n / 8);
        for (size_t i = 6; i > 0; --i) {
            calc::array_add(tmp_array_, tmp_array_, b_tab_[i - 1]);
            calc::array_scale(tmp_array_, tmp_array_, i * h_n / (i + 1));
        }
        calc::array_add(tmp_array_, tmp_array_, acceleration0_);
        calc::array_advance(tmp_state_, tmp_array_, step_size * h_n);

        particles.read_from_scalar_array(tmp_state_);
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussDadau<TypeSystem>::integrate_at_end(ParticleSys &particles, Scalar step_size) {
        tmp_state_ = input_;
        for (size_t i = 7; i > 0; --i) {
            calc::array_advance(tmp_state_, b_tab_[i - 1], step_size / (i + 1));
        }
        calc::array_advance(tmp_state_, acceleration0_, step_size);
        particles.read_from_scalar_array(tmp_state_);
    }

    template <typename TypeSystem>
    void GaussDadau<TypeSystem>::update_B_table(const ScalarArray &acc0, const ScalarArray &acc, size_t stage) {
        calc::array_sub(tmp_array_, acc, acc0);
        calc::array_scale(tmp_array_, tmp_array_, RadauConsts::G_tab(stage, 0));
        for (size_t j = 0; j < stage; ++j) {
            calc::array_advance(tmp_array_, g_tab_[j], -RadauConsts::G_tab(stage, j + 1));
        }
        calc::array_sub(dg_array_, tmp_array_, g_tab_[stage]);

        swap(g_tab_[stage], tmp_array_);

        /*calc::array_sub(tmp_array_, acc, acc0);
        calc::array_div_scale(tmp_array_, tmp_array_, RadauConsts::RR_tab(stage, 0));

        for (size_t j = 0; j < stage; ++j) {
            calc::array_retreat(tmp_array_, g_tab_[j]);
            calc::array_div_scale(tmp_array_, tmp_array_, RadauConsts::RR_tab(stage, j + 1));
        }

        calc::array_sub(dg_array_, tmp_array_, g_tab_[stage]);
        swap(g_tab_[stage], tmp_array_);*/

        for (size_t i = 0; i <= stage; ++i) {
            calc::array_advance(b_tab_[i], dg_array_, RadauConsts::G2B_tab(stage, i));
        }
    }

    template <typename TypeSystem>
    template <typename ParticleSys>
    void GaussDadau<TypeSystem>::calc_B_table(ParticleSys &particles, Scalar step_size) {
        particles.write_to_scalar_array(input_);
        check_particle_size(input_.size());
        particles.evaluate_general_derivative(acceleration0_);

        for (size_t i = 0; i < final_point; ++i) {
            integrate_to(particles, step_size, i);
            particles.evaluate_general_derivative(acceleration_);
            update_B_table(acceleration0_, acceleration_, i);
            particles.read_from_scalar_array(input_);
        }
    }

    template <typename TypeSystem>
    void GaussDadau<TypeSystem>::predict_new_B(Scalar step_ratio) {
        std::array<Scalar, final_point> Q;
        Q[0] = step_ratio;
        Q[1] = Q[0] * Q[0];
        Q[2] = Q[1] * Q[0];
        Q[3] = Q[1] * Q[1];
        Q[4] = Q[2] * Q[1];
        Q[5] = Q[2] * Q[2];
        Q[6] = Q[3] * Q[2];

        for (size_t i = 0; i < final_point; ++i) {
            calc::array_sub(tmp_array_, b_tab_[i], old_b_tab_[i]);
            calc::array_scale(old_b_tab_[i], b_tab_[6], Q[i] * RadauConsts::B_tab(6, i));
            for (size_t j = 6; j > i; --j) {
                calc::array_advance(old_b_tab_[i], b_tab_[j - 1], Q[i] * RadauConsts::B_tab(j - 1, i));
            }
            //  calc::array_scale(old_b_tab_[i], old_b_tab_[i], Q[i]);
            calc::array_add(b_tab_[i], old_b_tab_[i], tmp_array_);
        }

        RadauConsts::transfer_B_to_G(b_tab_, g_tab_);
    }
}  // namespace space::integrator
