#ifndef DYNAMICSYSTEM_H
#define DYNAMICSYSTEM_H

#include <functional>
#include "core-computation.hpp"
#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class RunArgs Declaration
\*---------------------------------------------------------------------------*/
  template<typename ParticleSys>
  class RunArgs {
  public:
    // type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(ParticleSys);

    using Callback = std::function<void(ParticleSys &)>;

    using StopCall = std::function<bool(ParticleSys &)>;

    // public members
    Scalar step_size{0};

    Scalar end_time{0};

    Scalar atol{6e-14};

    Scalar rtol{1e-13};

    // public methods
    /**
     *
     * @param partc_sys
     */
    void pre_operations(ParticleSys &partc_sys) const;

    /**
     *
     * @param partc_sys
     */
    void post_operations(ParticleSys &partc_sys) const;

    /**
     *
     * @param partc_sys
     */
    void stop_operations(ParticleSys &partc_sys) const;

    /**
     *
     * @param partc_sys
     * @return
     */
    bool check_stops(ParticleSys &partc_sys) const;

    /**
     *
     * @tparam Func
     * @tparam Args
     * @param func
     * @param args
     */
    template<typename Func, typename... Args>
    void add_pre_step_operation(Func func, Args &&... args);

    /**
     *
     * @tparam Func
     * @tparam Args
     * @param func
     * @param args
     */
    template<typename Func, typename... Args>
    void add_post_step_operation(Func func, Args &&... args);

    /**
     *
     * @tparam Func
     * @tparam Args
     * @param func
     * @param args
     */
    template<typename Func, typename... Args>
    void add_stop_point_operation(Func func, Args &&... args);

    /**
     *
     * @tparam Func
     * @tparam Args
     * @param func
     * @param args
     */
    template<typename Func, typename... Args>
    void add_stop_condition(Func func, Args &&... args);

    /**
     *
     * @tparam Scalar
     * @param end_
     */
    void add_stop_condition(Scalar end_);

    /**
     * @brief
     *
     * @return true
     * @return false
     */
    bool is_end_time_set() const { return is_end_time_set_; }

    bool is_stop_condition_set() const { return stop_cond_.size() > 0; }

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
     * @param t
     * @param partc
     */
    template<typename STL>
    Simulator(Scalar t, STL const &partc);

    /**
     *
     * @tparam T
     * @param t
     * @param p
     */
    template<typename... T>
    explicit Simulator(Scalar t, T const &... p);

    // Public methods
    /**
     *
     * @param arg
     */
    void run(RunArgs const &arg);

    virtual ~Simulator() = default; /**< @brief Default destructor, virtualize for inherent class*/
  private:
    // Private methods
    void advance_one_step();

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
  void RunArgs<ParticleSys>::pre_operations(ParticleSys &partc_sys) const {
    for (auto const &opt : pre_opts_) {
      opt(partc_sys);
    }
  }

  template<typename ParticleSys>
  void RunArgs<ParticleSys>::post_operations(ParticleSys &partc_sys) const {
    for (auto const &opt : post_opts_) {
      opt(partc_sys);
    }
  }

  template<typename ParticleSys>
  void RunArgs<ParticleSys>::stop_operations(ParticleSys &partc_sys) const {
    for (auto const &opt : stop_opts_) {
      opt(partc_sys);
    }
  }

  template<typename ParticleSys>
  bool RunArgs<ParticleSys>::check_stops(ParticleSys &partc_sys) const {
    for (auto const &check : stop_cond_) {
      if (check(partc_sys)) return true;
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
  void RunArgs<ParticleSys>::add_stop_condition(Scalar end_) {
    end_time = end_;
    is_end_time_set_ = true;
  }

/*---------------------------------------------------------------------------*\
    Class Simulator Implememtation
\*---------------------------------------------------------------------------*/
  template<typename ParticSys, typename OdeIterator>
  template<typename STL>
  Simulator<ParticSys, OdeIterator>::Simulator(Scalar t, const STL &partc) : particles_(t, partc) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
  }

  template<typename ParticSys, typename OdeIterator>
  template<typename... T>
  Simulator<ParticSys, OdeIterator>::Simulator(Scalar t, T const &... p)
          : Simulator(t, std::initializer_list<Particle>{p...}) {
    static_assert(calc::all(std::is_same_v<T, Particle>...), "Wrong particles type!");
  }

  template<typename ParticSys, typename OdeIterator>
  void Simulator<ParticSys, OdeIterator>::run(RunArgs const &arg) {
    if (!arg.is_stop_condition_set() && !arg.is_end_time_set()) {
      space::spacehub_abort("Use 'add_stop_condition' to set stop condition.");
    }

    step_size_ = arg.step_size;

    if (iseq(step_size_, 0.0)) {
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

    for (; particles_.time() < end_time;) {
      if (arg.check_stops(particles_)) break;

      arg.pre_operations(particles_);
      advance_one_step();
      arg.post_operations(particles_);
    }
    arg.stop_operations(particles_);
  }

  template<typename ParticSys, typename OdeIterator>
  void Simulator<ParticSys, OdeIterator>::advance_one_step() {
    particles_.pre_iter_process();
    step_size_ = iterator_.iterate(particles_, step_size_);
    particles_.post_iter_process();
  }

}  // namespace space
#endif
