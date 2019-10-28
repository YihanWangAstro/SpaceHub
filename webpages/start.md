@page start Getting started

@tableofcontents

@section prerequisite Prerequisite 
SpaceHub is written based on the morden-ish C++ standard @blabel{c++17}. To make the code compiled, the compiler on the target
machine need to give full support on at least @blabel{c++17}. The lowest compiler version for different OS is listed below.

| Operating System | Compiler             | Lowest Version required | Version recommended    |
| ---------------- | -------------------- | ----------------------- | ---------------------- |
| Linux            | GCC                  | @bflabel{7}             | @gflabel{9}            |
| Linux/MacOS      | Intel C++            | @bflabel{19.01}         | @gflabel{19.01}        |
| MacOS            | Clang                | @bflabel{4}             | @gflabel{8}            |
| MacOS            | Apple Clang          | @bflabel{support}       | @gflabel{support}      |
| Windows          | MSVC                 | @bflabel{19.14}         | @gflabel{19.14}        |
| Windows          | GCC via CygWin/MinGw | @bflabel{7}             | @gflabel{9 if support} |

<ul>
<li>
For Linux/MacOS user, check the version of your GCC in terminal,

@code{.shell-session}
> g++ -v
@endcode

or check the version of the intel++ via,

@code{.shell-session}
> icpc -v
@endcode

If your compiler version doesn't satisfy the lowest version required by SpaceHub, you need to update/install the new compiler.

You can install the new GCC from scratch by executing the following commands by sequence,

@code{.shell-session}
> wget https://ftp.gnu.org/gnu/gcc/gcc-9.1.0/gcc-9.1.0.tar.xz
> tar xf gcc-9.1.0.tar.xz
> cd gcc-9.1.0
> ./contrib/download_prerequisites
> mkdir build
> cd build
> ../configure --prefix=$HOME/GCC-9.1.0 --enable-languages=c,c++ --disable-multilib
> make -j 8
> make install
@endcode

after the installation, add the new g++ executable to your environment variable `PATH` by adding the following lines into your `.bashrc`(Linux)/`.bash_profile`(MacOS) under your directory `~`,

@code
export PATH=$PATH:$HOME/GCC-9.1.0/bin
export LD_LIBRARY_PATH=$HOME/GCC-9.1.0/lib
export LD_LIBRARY_PATH=$HOME/GCC-9.1.0/lib64
@endcode

and source it

@code{.shell-session}
> source ~/.bashrc
@endcode

or

@code{.shell-session}
> source ~/.bash_profile
@endcode

then, check the GCC version again via

@code{.shell-session}
> g++ -v
@endcode

</li>

<li>
For Windows user, if you have installed [Microsoft Visual Studio](https://visualstudio.microsoft.com/vs/), check the version of your MSVC in cmd terminal,

@code{.shell-session}
path_to_your_installed_dir\Microsoft Visual Studio xxx\VC > cl
@endcode

or just open your [Microsoft Visual Studio](https://visualstudio.microsoft.com/vs/), check the version information through _about_ or _help_.

If your MSVC doesn't satisfy the lowest version required by SpaceHub, update [Microsoft Visual Studio](https://visualstudio.microsoft.com/vs/) to the newest version.
If [Microsoft Visual Studio](https://visualstudio.microsoft.com/vs/) is too heavy for you, you can also install the GCC on Windows via [CygWin](http://www.cygwin.com/) or [MinGw](http://www.mingw.org/).
</li>

<li>
For cluster user, if the cluster has build system like [EasyBuild](https://easybuild.readthedocs.io/en/latest/), check if any gcc module with version later than 7 is available via


@code{.shell-session}
> module avail
@endcode

If any GCC module is available, e.g @blabel{gcc/9.2.0}, you can use 

@code{.shell-session}
> module load gcc/9.2.0
@endcode

to load the GCC-9.2.0. Otherwise, you need to install the GCC under your own directories without root privilege through the method above for linux user. 

</li>
</ul>

@section install Install 
SpaceHub is a header only code(In @blabel{c++20}, we will make into [Modules](https://en.cppreference.com/w/cpp/language/modules)). @m_span{m-text m-danger} No installation @m_endspan is required, you only need
@code{.cpp}#include"path_to_spacehub/src/spaceHub.hpp" 
@endcode in your application code to use the SpaceHub.
<ul>
<li>
For Linux/MacOS user, open the terminal and `cd` to a place where you want to install the SpaceHub,


@code{.shell-session}
> git clone https://github.com/YihanWangAstro/SpaceHub.git
@endcode

That's it.
</li>

<li>
</li>
For Windows user, download the repository of SpaceHub and unpack it.
@m_div{m-button m-info} <a href="https://github.com/YihanWangAstro/SpaceHub/archive/master.zip">@m_div{m-big}SpaceHub @m_enddiv @m_div{m-small} download @m_enddiv </a> @m_enddiv 
</ul>

@section demo Demo 
Here is a very simple example to integrate a (sun-earth-moon) system. For more information we recommend to read tutorial carefully.

@code{.cpp}
// main.cpp
#include"PATH_TO_SPACEHUB/src/spaceHub.hpp"
using namespace space;
using namespace unit;
int main(){
  using Particle = typename DefaultSolver::Particle;

  // create three particles. particles are rest at original point.
  Particle sun{m_solar}, earth{m_earth}, moon{m_moon};

  // create a Kepler orbit of (moon mass, earth mass) with a = 268782 km, e = 0.055 and i = 1.543 degree.
  auto moon_orbit = orbit::EllipOrbit{earth.mass, moon.mass, 268782 * km, 0.055, 1.543 * deg, 0.0, 0.0, 0.0};

  // move the moon to the moon orbit
  orbit::move_particles_to(moon_orbit, moon);

  // create a Kepler orbit of (moon mass + earth mass, solar mass) a = 1 au, e = 0.016 and i = 7.155 degree.
  auto earth_orbit = orbit::EllipOrbit{sun.mass, earth.mass + moon.mass, au, 0.016, 7.155 * deg, 0.0, 0.0, 0.0};

  // move the centre of mass of the moon and earth to the earth orbit.
  orbit::move_particles_to(earth_orbit, earth, moon);

  // move the three objects to the centre of mass reference frame
  orbit::move_to_com_coord(sun, earth, moon);

  // Initialize the system with the three particles and set the time = 0.
  DefaultSolver sim{0, sun, earth, moon};

  DefaultSolver::RunArgs args;

  // add a default output printer
  args.add_pre_step_operation(run_operations::DefaultWriter("solar.dat"));

  args.add_stop_condition(100 * unit::year);

  // run simulation with arguments.
  sim.run(args);

  return 0;
}
@endcode

compile it with

@code{.ansi}
> g++ -std=c++17 -O3 -o test main.cpp
@endcode

and run 


@code{.ansi}
> ./test
@endcode

Enjoy!

