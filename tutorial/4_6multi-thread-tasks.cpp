
#include "../src/spaceHub.hpp"
#include "../src/taskflow/taskflow.hpp"

using namespace space;
using namespace unit;
using namespace callback;
using namespace orbit;
using Solver = methods::DefaultMethod<>;
using Particle = Solver::Particle;

void job(size_t task_id, size_t scattering_num) {
    double v_inf = 10_kms;
    double b_max = 10_AU;
    double r_start = 100_AU;

    std::fstream file("scattering_par_" + std::to_string(task_id) + ".txt", std::ios::out);

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle p1{1_Ms};
        Particle p2{1_Ms};

        auto orb = scattering::incident_orbit(p1.mass, p2.mass, v_inf, b_max, r_start);

        orbit::move_particles(orb, p2);

        orbit::move_to_COM_frame(p1, p2);

        Solver solver{0, p1, p2};

        Solver::RunArgs args;

        double t_end = 2 * orbit::time_to_periapsis(p1, p2);

        args.add_stop_condition(t_end);

        args.add_stop_point_operation([&file, &orb](auto& ptc, auto h) { print(file, orb, ',', ptc, '\n'); });

        solver.run(args);
    }
    print(std::cout, "job ", task_id, " finished\n");
}

int main(int argc, char** argv) {
    size_t n = 500000;  // total scattering number
    // job number; any number > cpu core number can utilize all CPU powers; jobs will be scheducled by executor so don't
    // worry about the cpu resource competition
    size_t job_num = 24;

    tf::Executor executor;  // create multi thread executor; thread number = cpu logical core number

    tools::Timer timer;
    timer.start();
    for (size_t i = 0; i < job_num; ++i) {
        // put jobs in job queue. They will be scheduled and executed by executor in parallel
        executor.silent_async(job, i, n / job_num);
    }
    executor.wait_for_all();  // wait all jobs to be finished
    print(std::cout, "finished with ", multi_thread::machine_thread_num, " threads in : ", timer.get_time(), " s\n");

    return 0;
}
