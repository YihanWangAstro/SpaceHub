#ifndef SPACEHUB_SIMULATOR_HPP
#define SPACEHUB_SIMULATOR_HPP

#include <functional>
#include "core-computation.hpp"
#include "dev-tools.hpp"


namespace space {

/*---------------------------------------------------------------------------*\
    Class RunArgs Declaration
\*---------------------------------------------------------------------------*/
/**
 * Run arguments that is used to set all arguments needed by Simulator.
 *
 * @tparam ParticleSys Any implementation of concept ParticleSystem.
 */
  template<typename ParticleSys>
  class RunArgs {
  public:
    // type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);

    /**
     * Callback function type for pre-operation, pos-operation and stop operation.
     */
    using Callback = std::function<void(ParticleSys &)>;

    /**
     * Callback function type for stop condition.
     */
    using StopCall = std::function<bool(ParticleSys &)>;

    // public members

    /**
     * Initial step size for the integration of the Simulator.
     */
    Scalar step_size{0};

    /**
     * The time duration for the integration of the Simulator.
     */
    Scalar end_time{0};

    /**
     * The absolute error tolerance.
     */
    Scalar atol{1e-12};

    /**
     * The relative error tolerance.
     */
    Scalar rtol{1e-12};

    // public methods
    /**
     * Call the all registered pre-operation functions by sequence.
     *
     * @param[in] particle_system The particle system that is going to be operated.
     */
    void pre_operations(ParticleSys &particle_system) const;

    /**
     * Call the all registered pos-operation functions by sequence.
     *
     * @param[in] particle_system The particle system that is going to be operated.
     */
    void post_operations(ParticleSys &particle_system) const;

    /**
     * Call the all registered stop-operation functions by sequence.
     *
     * @param[in] particle_system The particle system that is going to be operated.
     */
    void stop_operations(ParticleSys &particle_system) const;

    /**
     * Check the all registered stop condition functions by sequence. If any of them is satisfied, return `true`.
     *
     * @param[in] particle_system The particle system that is going to be checked
     */
    bool check_stops(ParticleSys &particle_system) const;

    /**
     * Register a callable object(function pointer, functor, lambda,etc...) to pre-step-operations.
     *
     * @tparam Func Callable type that is conversional to `Callback`.
     * @tparam Args Type of the binding arguments.
     * @param[in] func Callable object.
     * @param[in] args Binding arguments.
     */
    template<typename Func, typename... Args>
    void add_pre_step_operation(Func func, Args &&... args);

    /**
     * Register a callable object(function pointer, functor, lambda,etc...) to post-step-operations.
     *
     * @tparam Func Callable type that is conversional to `Callback`.
     * @tparam Args Type of the binding arguments.
     * @param[in] func Callable object.
     * @param[in] args Binding arguments.
     */
    template<typename Func, typename... Args>
    void add_post_step_operation(Func func, Args &&... args);

    /**
     * Register a callable object(function pointer, functor, lambda,etc...) to stop-point-operations.
     *
     * @tparam Func Callable type that is conversional to `Callback`.
     * @tparam Args Type of the binding arguments.
     * @param[in] func Callable object.
     * @param[in] args Binding arguments.
     */
    template<typename Func, typename... Args>
    void add_stop_point_operation(Func func, Args &&... args);

    /**
     * Register a callable object(function pointer, functor, lambda,etc...) to stop conditions.
     *
     * @tparam Func Callable type that is conversional to `StopCall`.
     * @tparam Args Type of the binding arguments.
     * @param[in] func Callable object.
     * @param[in] args Binding arguments.
     */
    template<typename Func, typename... Args>
    void add_stop_condition(Func func, Args &&... args);

    /**
     * Add the duration time of the integration as a stop condition.
     * @tparam T Floating point like scalar.
     * @param[in] end Duration time of the integration.
     */
    template<typename T>
    void add_stop_condition(T end);

    /**
     * Check if the integration duration time is set.
     * @return boolean
     */
    [[nodiscard]] bool is_end_time_set() const { return is_end_time_set_; }

    /**
     * Check if any of the stop condition(except the duration time) is set.
     * @return boolean
     */
    [[nodiscard]] bool is_stop_condition_set() const { return stop_cond_.size() > 0; }

  private:
    // private members
    std::vector<Callback> pre_opts_;

    std::vector<Callback> post_opts_;

    std::vector<Callback> stop_opts_;

    std::vector<StopCall> stop_cond_;

    bool is_end_time_set_{false};
  };

/*---------------------------------------------------------------------------*\
    Class Simulator Declaration
\*---------------------------------------------------------------------------*/
/**
 * 
 * @tparam ParticSys
 * @tparam OdeIterator
 */
  template<typename ParticSys, typename OdeIterator>
  class Simulator {
  public:
    // Type member
    SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);

    using RunArgs = space::RunArgs<ParticSys>;

    using Particle = typename ParticSys::Particle;

    SPACEHUB_READ_ACCESSOR(auto, particles, particles_);

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(Simulator, delete, default, default, default, default);

    /**
     *
     * @tparam STL
     * @param time
     * @param particle_set
     */
    template<typename STL>
    Simulator(Scalar time, STL const &particle_set);

    /**
     *
     * @tparam T
     * @param time
     * @param particle
     */
    template<typename... T>
    explicit Simulator(Scalar time, T const &... particle);

    // Public methods
    /**
     *
     * @param arg
     */
    void run(RunArgs const &arg);

    virtual ~Simulator() = default; /**< @brief Default destructor, virtualize for inherent class*/
  private:
    // Private methods
    inline void advance_one_step();

    // Private members
    /** @brief Macro step size for ODE iterator*/
    Scalar step_size_{0.0};

    /** @brief Particle system*/
    ParticSys particles_;

    /** @brief ODE Iterator*/
    OdeIterator iterator_;

    CREATE_CRTP_IMPLEMENTATION_CHECK(set_atol);

    CREATE_CRTP_IMPLEMENTATION_CHECK(set_rtol);
  };

/*---------------------------------------------------------------------------*\
    Class RunArgs Implementation
\*---------------------------------------------------------------------------*/
  template<typename ParticleSys>
  void RunArgs<ParticleSys>::pre_operations(ParticleSys &particle_system) const {
    for (auto const &opt : pre_opts_) {
      opt(particle_system);
    }
  }

  template<typename ParticleSys>
  void RunArgs<ParticleSys>::post_operations(ParticleSys &particle_system) const {
    for (auto const &opt : post_opts_) {
      opt(particle_system);
    }
  }

  template<typename ParticleSys>
  void RunArgs<ParticleSys>::stop_operations(ParticleSys &particle_system) const {
    for (auto const &opt : stop_opts_) {
      opt(particle_system);
    }
  }

  template<typename ParticleSys>
  bool RunArgs<ParticleSys>::check_stops(ParticleSys &particle_system) const {
    for (auto const &check : stop_cond_) {
      if (check(particle_system)) return true;
    }
    return false;
  }

  template<typename ParticleSys>
  template<typename Func, typename... Args>
  void RunArgs<ParticleSys>::add_pre_step_operation(Func func, Args &&... args) {
    pre_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
  }

  template<typename ParticleSys>
  template<typename Func, typename... Args>
  void RunArgs<ParticleSys>::add_post_step_operation(Func func, Args &&... args) {
    post_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
  }

  template<typename ParticleSys>
  template<typename Func, typename... Args>
  void RunArgs<ParticleSys>::add_stop_point_operation(Func func, Args &&... args) {
    stop_opts_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
  }

  template<typename ParticleSys>
  template<typename Func, typename... Args>
  void RunArgs<ParticleSys>::add_stop_condition(Func func, Args &&... args) {
    stop_cond_.emplace_back(std::bind(std::forward<Func>(func), std::placeholders::_1, std::forward<Args>(args)...));
  }

  template<typename ParticleSys>
  template<typename T>
  void RunArgs<ParticleSys>::add_stop_condition(T end) {
    end_time = end;
    is_end_time_set_ = true;
  }

/*---------------------------------------------------------------------------*\
    Class Simulator Implememtation
\*---------------------------------------------------------------------------*/
  template<typename ParticSys, typename OdeIterator>
  template<typename STL>
  Simulator<ParticSys, OdeIterator>::Simulator(Scalar time, const STL &particle_set) : particles_(time, particle_set) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
  }

  template<typename ParticSys, typename OdeIterator>
  template<typename... T>
  Simulator<ParticSys, OdeIterator>::Simulator(Scalar time, T const &... particle)
          : Simulator(time, std::initializer_list<Particle>{particle...}) {
    static_assert(calc::all(std::is_same_v<T, Particle>...), "Wrong particles type!");
  }

  template<typename ParticSys, typename OdeIterator>
  void Simulator<ParticSys, OdeIterator>::run(RunArgs const &arg) {
    if (!arg.is_stop_condition_set() && !arg.is_end_time_set()) {
      space::spacehub_abort("Use 'add_stop_condition' to set stop condition.");
    }

    step_size_ = arg.step_size;

    if (step_size_ == 0.0) {
      step_size_ = 0.01 * calc::calc_step_scale(particles_);
    }

    Scalar end_time = space::unit::hubble_t;

    if (arg.is_end_time_set()) {
      end_time = arg.end_time;
    }

    if (particles_.time() >= end_time) {
      space::print(std::cout, "Warning: The stop time is '<=' to the start time!");
    }

    if constexpr (HAS_CRTP_IMPLEMENTATION(OdeIterator, set_atol, Scalar)) {
      iterator_.set_atol(arg.atol);
    }

    if constexpr (HAS_CRTP_IMPLEMENTATION(OdeIterator, set_rtol, Scalar)) {
      iterator_.set_rtol(arg.rtol);
    }

    for (; particles_.time() < end_time && !arg.check_stops(particles_);) {

      arg.pre_operations(particles_);
      advance_one_step();
      arg.post_operations(particles_);
    }
    arg.stop_operations(particles_);
  }

  template<typename ParticSys, typename OdeIterator>
  inline void Simulator<ParticSys, OdeIterator>::advance_one_step() {
    particles_.pre_iter_process();
    step_size_ = iterator_.iterate(particles_, step_size_);
    particles_.post_iter_process();
  }

}  // namespace space
#endif
