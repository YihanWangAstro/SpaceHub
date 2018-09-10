#ifndef TEST_TOOLS_H
#define TEST_TOOLS_H
#include "type_class.h"
namespace UnitTest {
    template <size_t ArraySize, typename T>
    void print_type_base() {
        if constexpr (ArraySize == SpaceH::DYNAMICAL)
            std::cout << "Dynamical Array";
        else
            std::cout << "Fixed size Array";

        std::cout << " " << typeid(T).name() << '\n';
    }
}
#endif