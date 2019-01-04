
#ifndef PARTICLES_H
#define PARTICLES_H

#include "dev_tools.h"

namespace SpaceH {


    template<typename TypeSystem>
    class BasicParticles {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(px, ScalarArray, px_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(py, ScalarArray, py_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(pz, ScalarArray, pz_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vx, ScalarArray, vx_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vy, ScalarArray, vy_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vz, ScalarArray, vz_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(mass, ScalarArray, mass_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(idn, IndexArray, idn_);

        BasicParticles() = default;

        template<typename Container>
        explicit BasicParticles(Container const &partc) {
            size_t input_num = partc.size();
            if (input_num > _capacity_) {
                SPACEHUB_ABORT("Input particle number exceeds the capacity!");
            } else {
                size_t i = 0;
                for(auto& p : partc){
                    std::tie(idn(i), mass(i), px(i), py(i), pz(i), vx(i), vy(i), vz(i)) = p.basic_info();
                    i++;
                }
                active_num = input_num;
            }
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, BasicParticles const &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::printd(' ', os, ps.idn(i), ps.mass(i), ps.px(i), ps.py(i), ps.pz(i), ps.vx(i), ps.vy(i), ps.vz(i));
            }
        }

        friend std::istream &operator>>(std::istream &is, BasicParticles &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn(i), ps.mass(i), ps.px(i), ps.py(i), ps.pz(i), ps.vx(i), ps.vy(i), ps.vz(i));
            }
        }

    private:
        ScalarArray px_;
        ScalarArray py_;
        ScalarArray pz_;
        ScalarArray vx_;
        ScalarArray vy_;
        ScalarArray vz_;
        ScalarArray mass_;
        IndexArray idn_;
        size_t active_num{_capacity_};
    };


}
#endif

