
#ifndef GAUSSDADAU_H
#define GAUSSDADAU_H
#include "src/dev_tools.h"
#include "src/type_class.h"
namespace SpaceH
{
    namespace Radau
    {
        /** Gauss radau stepping. see details in
         https://www.cambridge.org/core/journals/international-astronomical-union-colloquium/article/an-efficient-integrator-that-uses-gauss-radau-spacings/F942BC9121C74CC2FA296050FC18D824
         */

        constexpr double interplt[8] = {5.62625605369221465e-02, 1.80240691736892365e-01, 3.52624717113169637e-01,   5.47153626330555383e-01, 7.34210177215410532e-01, 8.85320946839095768e-01, 9.77520613561287502e-01};

        constexpr double VC[8][8] =
        {
            {5.62625605369221488e-02, 1.58273785908541475e-03, 5.93659230739144700e-05, 2.50505913058228221e-06, 1.12752832786364175e-07, 5.28646923362689570e-09, 2.54940253100151135e-10, 1.25506424954273193e-11},
            {1.80240691736892361e-01, 1.62433534788967291e-02, 1.95180884477546898e-03, 2.63846532240386473e-04, 3.80447051867100241e-05, 5.71433664981562616e-06, 8.82819420497352435e-07, 1.39229985150554597e-07},
            {3.52624717113169617e-01, 6.21720955595714486e-02, 1.46156117393512206e-02, 3.86536946626848384e-03, 1.09041985166464603e-03, 3.20424159773191819e-04, 9.68481245967829960e-05, 2.98821622215214072e-05},
            {5.47153626330555420e-01, 1.49688545403338535e-01, 5.46017536250551146e-02, 2.24066606249673420e-02, 9.80790849192716036e-03, 4.47202724839682903e-03, 2.09733079372232617e-03, 1.00411688072492376e-03},
            {7.34210177215410487e-01, 2.69532292163342236e-01, 1.31928901329682199e-01, 7.26476565188252785e-02, 4.26709190135767829e-02, 2.61078525090855322e-02, 1.64302723006367039e-02, 1.05553639953544340e-02},
            {8.85320946839095790e-01, 3.91896589456036537e-01, 2.31302839760153781e-01, 1.53582936827273248e-01, 1.08776152840200469e-01, 8.02515055216670650e-02, 6.08985761603187514e-02, 4.71754369689804106e-02},
            {9.77520613561287499e-01, 4.77773274968617985e-01, 3.11355483260339450e-01, 2.28267302274238649e-01, 1.78508794700074913e-01, 1.45413355434419279e-01, 1.21838187792222097e-01, 1.04211922575117274e-01},
            {1.00000000000000000e+00, 5.00000000000000000e-01, 3.33333333333333333e-01, 2.50000000000000000e-01, 2.00000000000000000e-01, 1.66666666666666667e-01, 1.42857142857142857e-01, 1.25000000000000000e-01}
        };

        constexpr double PC[8][9] =
        {
            {5.62625605369221488e-02, 1.58273785908541475e-03, 2.96829615369572350e-05, 8.35019710194094070e-07, 2.81882081965910438e-08, 1.05729384672537914e-09, 4.24900421833585225e-11, 1.79294892791818847e-12, 7.84590314640274681e-14},
            {1.80240691736892361e-01, 1.62433534788967291e-02, 9.75904422387734491e-04, 8.79488440801288243e-05, 9.51117629667750602e-06, 1.14286732996312523e-06, 1.47136570082892073e-07, 1.98899978786506567e-08, 2.78832320378369023e-09},
            {3.52624717113169617e-01, 6.21720955595714486e-02, 7.30780586967561031e-03, 1.28845648875616128e-03, 2.72604962916161508e-04, 6.40848319546383638e-05, 1.61413540994638327e-05, 4.26888031736020103e-06, 1.17079877778820337e-06},
            {5.47153626330555420e-01, 1.49688545403338535e-01, 2.73008768125275573e-02, 7.46888687498911401e-03, 2.45197712298179009e-03, 8.94405449679365805e-04, 3.49555132287054362e-04, 1.43445268674989108e-04, 6.10451325053742022e-05},
            {7.34210177215410487e-01, 2.69532292163342236e-01, 6.59644506648410994e-02, 2.42158855062750928e-02, 1.06677297533941957e-02, 5.22157050181710644e-03, 2.73837871677278399e-03, 1.50790914219349057e-03, 8.61095074400260249e-04},
            {8.85320946839095790e-01, 3.91896589456036537e-01, 1.15651419880076890e-01, 5.11943122757577493e-02, 2.71940382100501173e-02, 1.60503011043334130e-02, 1.01497626933864586e-02, 6.73934813842577294e-03, 4.64060028054731337e-03},
            {9.77520613561287499e-01, 4.77773274968617985e-01, 1.55677741630169725e-01, 7.60891007580795496e-02, 4.46271986750187282e-02, 2.90826710868838558e-02, 2.03063646320370162e-02, 1.48874175107310391e-02, 1.13188113884477806e-02},
            {1.00000000000000000e+00, 5.00000000000000000e-01, 1.66666666666666667e-01, 8.33333333333333333e-02, 5.00000000000000000e-02, 3.33333333333333333e-02, 2.38095238095238095e-02, 1.78571428571428571e-02, 1.38888888888888889e-02}
        };

        /*constexpr double r[28] =
         {
             1.77738089140780001e+01,
             5.54813671853721662e+00, 8.06593864838188709e+00,
             2.83587607864443884e+00, 3.37424997696263552e+00, 5.80100155926406205e+00,
             1.82764026751759770e+00, 2.03711183535858464e+00, 2.72544221180822598e+00, 5.14062410581093270e+00,
             1.36200781606246958e+00, 1.47504021756041165e+00, 1.80515358014025139e+00, 2.62064492638703525e+00, 5.34597689987110987e+00,
             1.12953387533678987e+00, 1.20618766605844559e+00, 1.41827826373473910e+00, 1.87724249618680995e+00, 2.95711601729045588e+00, 6.61766201370242155e+00,
             1.02299632982348675e+00, 1.08547219393864239e+00, 1.25426462228187776e+00, 1.60026654949081622e+00, 2.32359830021969444e+00, 4.10997577834455837e+00, 1.08460261902368476e+01
         };*/

        /*coef for computing g(derived from r)*/
        constexpr double gg[28] =
        {
            1.77738089140780001e+01,
            4.47509303845559954e+01, 8.06593864838188709e+00,
            5.55095216749226992e+01, 1.95740293777069741e+01, 5.80100155926406205e+00,
            5.21625022561530086e+01, 2.85409022679298973e+01, 1.40104739330160323e+01, 5.14062410581093270e+00,
            5.08080910907447823e+01, 3.73038175637124453e+01, 2.52900342103279856e+01, 1.40099072392295156e+01, 5.34597689987110987e+00,
            7.09853803416487244e+01, 6.28448441358017857e+01, 5.21020450666394487e+01, 3.67361232269326421e+01, 1.95691943377340431e+01, 6.61766201370242155e+00,
            2.30858165231426693e+02, 2.25668615322657310e+02, 2.07899029180855735e+02, 1.65753721732680304e+02, 1.03578820531755138e+02, 4.45769049331641530e+01, 1.08460261902368476e+01
        };

        constexpr double c[21] =
        {
             1.27179030902686775e-03, -1.43653023637089149e-03,  1.95656540994722109e-03, -3.57589772925161718e-03,  1.01408028300636298e-02, -5.62625605369221488e-02,
            -3.87603579159067693e-02,  4.21585277212687057e-02, -5.47553868890686898e-02,  9.35376952594620670e-02, -2.36503252273814524e-01,
             3.60962243452845999e-01, -3.60099596502056807e-01,  4.15881200082306890e-01, -5.89127969386984196e-01,
            -1.46688420840042699e+00,  1.25015071184069093e+00, -1.13628159571753962e+00,
             2.90613625930842900e+00, -1.87049177293294999e+00,
            -2.75581271977204567e+00
        };

        template<typename RadauTab>
        void transG2B(const RadauTab& G, RadauTab& B)
        {
            size_t size = G.size();
            for(size_t i = 0 ; i < size; ++i) {
                B[i][0] = c[0] *G[i][6] + c[1] *G[i][5] + c[2] *G[i][4] + c[3] *G[i][3] + c[4] *G[i][2] + c[5]*G[i][1] + G[i][0];
                B[i][1] = c[6] *G[i][6] + c[7] *G[i][5] + c[8] *G[i][4] + c[9] *G[i][3] + c[10]*G[i][2] +      G[i][1];
                B[i][2] = c[11]*G[i][6] + c[12]*G[i][5] + c[13]*G[i][4] + c[14]*G[i][3] +       G[i][2];
                B[i][3] = c[15]*G[i][6] + c[16]*G[i][5] + c[17]*G[i][4] +       G[i][3];
                B[i][4] = c[18]*G[i][6] + c[19]*G[i][5] +       G[i][4];
                B[i][5] = c[20]*G[i][6] +       G[i][5];
                B[i][6] =       G[i][6];
            }
        }
    }

    /** @brief Gauss-Dadau integrator */
    template <typename ParticSys>
    class GaussDadau
    {
    public:
        /* Typedef */
        using type        = typename ParticSys::type;
        using Scalar      = typename type::Scalar;
        using Vector      = typename type::Vector;
        using VectorArray = typename type::VectorArray;

        template<typename T, size_t S>
        using Container = typename type::template Container<T, S>;

        using RadauArray = Container<Vector,7>;
        using RadauTab   = Container<RadauArray, ParticSys::arraySize>;
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        /** @brief Order of the integrator*/
        static const int order{15};
        static const size_t finalPoint{7};

        /** @brief Interface to integrate particle system
         *
         *  This function integrate the particle system for one step with Gauss-Radau stepping.
         *  @param particles  Particle system need to be integrated.
         *  @param stepLength Step size for integration.
         */
        void integrate(ParticSys& particles, Scalar stepLength) {
            checkTabVolume(particles.particleNumber());
            calcuBTab(particles, stepLength);
            evaluateSystemAt(particles, stepLength, finalPoint);
        }

        /**
         *
         * @param particles
         * @param stepLength
         */
        void calcuBTab(const ParticSys& particles, Scalar stepLength) {
            for(size_t i = 0 ; i < finalPoint; ++i) {
                localSystem_ = particles;
                evaluateSystemAt(localSystem_, stepLength, i);
                calcuGTab(particles.acc(), localSystem_.acc(), i);
            }
            Radau::transG2B(Gtab_, Btab_);
        }

        /**
         *
         * @param particleSys
         * @param stepLength
         * @param index
         */
        void evaluateSystemAt(ParticSys& particleSys, Scalar stepLength, size_t index) {
            VectorArray dpos;
            VectorArray dvel;
            evaluateVelIncrement(dvel, particleSys.acc(), index);
            evaluatePosIncrement(dpos, particleSys.vel(), particleSys.acc(), stepLength, index);
            particleSys.advanceVel(dvel, stepLength);
            particleSys.advancePos(dpos, stepLength);
            particleSys.advanceTime(stepLength);
            particleSys.evaluateAcc();
        }

        /**
         *
         * @return
         */
        inline const RadauTab& getBTab() const {
            return Btab_;
        }

        /**
         *
         * @return
         */
        inline const VectorArray& localAcc() const {
            return localSystem_.acc();
        }

        void checkTabVolume(size_t particleNum) {
            if constexpr (ParticSys::arraySize == SpaceH::DYNAMICAL){

                if (particleNum != particleNumber_) {

                    Btab_.resize(particleNum);
                    Gtab_.resize(particleNum);
                    oldBtab_.resize(particleNum);

                    //once particle number changes, old b Value should be droped away.
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

        /**
         *
         * @param Q1
         */
        void predictNewB(Scalar Q1) {
            Scalar Q2 = Q1*Q1;
            Scalar Q3 = Q2*Q1;
            Scalar Q4 = Q2*Q2;
            Scalar Q5 = Q3*Q2;
            Scalar Q6 = Q3*Q3;
            Scalar Q7 = Q4*Q3;
            size_t size = Btab_.size();
            RadauArray dB;
            for(size_t i = 0 ; i < size; ++i) {

                for(size_t j = 0 ; j < finalPoint; ++j)
                    dB[j] = Btab_[i][j] - oldBtab_[i][j];

                oldBtab_[i][0] = Q1*( 7*Btab_[i][6] +  6*Btab_[i][5] +  5*Btab_[i][4] +  4*Btab_[i][3] +  3*Btab_[i][2] +  2*Btab_[i][1] +  Btab_[i][0]);
                oldBtab_[i][1] = Q2*(21*Btab_[i][6] + 15*Btab_[i][5] + 10*Btab_[i][4] +  6*Btab_[i][3] +  3*Btab_[i][2] +    Btab_[i][1]);
                oldBtab_[i][2] = Q3*(35*Btab_[i][6] + 20*Btab_[i][5] + 10*Btab_[i][4] +  4*Btab_[i][3] +    Btab_[i][2]);
                oldBtab_[i][3] = Q4*(35*Btab_[i][6] + 15*Btab_[i][5] +  5*Btab_[i][4] +    Btab_[i][3]);
                oldBtab_[i][4] = Q5*(21*Btab_[i][6] +  6*Btab_[i][5] +    Btab_[i][4]);
                oldBtab_[i][5] = Q6*( 7*Btab_[i][6] +    Btab_[i][5]);
                oldBtab_[i][6] = Q7*    Btab_[i][6];

                for(size_t j = 0 ; j < 7; ++j)
                   Btab_[i][j] = oldBtab_[i][j] + dB[j];
            }
        }
    private:
        void evaluatePosIncrement(VectorArray& dpos, const VectorArray& vel, const VectorArray& acc, Scalar stepLength, size_t iter) {
            using namespace Radau;
            size_t size = vel.size();
            dpos.resize(size);

            for(size_t i = 0 ; i < size; ++i) {
                dpos[i] = vel[i]   *  PC[iter][0] + (acc[i]  *  PC[iter][1] + Btab_[i][0]*PC[iter][2] + Btab_[i][1]*PC[iter][3] + Btab_[i][2]*PC[iter][4]
                        + Btab_[i][3]*PC[iter][5] + Btab_[i][4]*PC[iter][6] + Btab_[i][5]*PC[iter][7] + Btab_[i][6]*PC[iter][8])*stepLength;
            }
        }

        void evaluateVelIncrement(VectorArray& dvel, const VectorArray& acc, size_t iter)
        {
            using namespace Radau;
            size_t size = acc.size();
            dvel.resize(size);
            for(size_t i = 0 ; i < size; ++i) {
                dvel[i] = acc[i]   *  VC[iter][0] + Btab_[i][0]*VC[iter][1] + Btab_[i][1]*VC[iter][2] + Btab_[i][2]*VC[iter][3]
                        + Btab_[i][3]*VC[iter][4] + Btab_[i][4]*VC[iter][5] + Btab_[i][5]*VC[iter][6] + Btab_[i][6]*VC[iter][7];
            }
        }

        /**
         *
         * @param a0
         * @param a
         * @param iter
         */
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
        RadauTab  oldBtab_;
        RadauTab  Btab_;
        RadauTab  Gtab_;
        ParticSys localSystem_;
        size_t    particleNumber_{ParticSys::arraySize};
    };

}
#endif
