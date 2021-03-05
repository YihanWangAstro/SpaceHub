
#include "../src/spaceHub.hpp"
using namespace space;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

int main(int argc, char **argv) {
    // mass, x , y, z, vx, vy, vz //units: AU = 1, year = 2pi, G = 1
    Particle sun{1.00000597682,        -4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6,
                 +6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9};
    Particle jup{1. / 1047.355,       +3.40546614227466e+0, +3.62978190075864e+0,   +3.42386261766577e-2,
                 -0.3254242234844626, 0.32078376079804843,  -0.00015504584274327233};
    Particle sat{1. / 3501.6,          +6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1,
                 -0.24261807906125546, 0.23236917361652312,  0.0009720111543216125};
    Particle ura{1. / 22869.,         +1.11636331405597e+1, +1.60373479057256e+1,  +3.61783279369958e-1,
                 -0.1894447922304644, 0.12000768830940597,  -0.0012655376696374542};
    Particle nep{1. / 19314.,          -3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1,
                 -0.01264216568441132, -0.18100181375010865, 0.0020831452402001265};

    orbit::move_to_COM_frame(sun, jup, sat, ura, nep);

    Solver solver{0, sun, jup, sat, ura, nep};  // first parameter t_start = 0;

    Solver::RunArgs args;

    args.add_stop_condition(1000 * unit::year);

    solver.run(args);

    print(std::cout, "first run without any output complete!\n");

    return 0;
}
