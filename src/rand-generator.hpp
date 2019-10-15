/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file rand-generator.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_RAND_GENERATOR_HPP
#define SPACEHUB_RAND_GENERATOR_HPP

#include <mutex>
#include <random>
#include "own-math.hpp"

/**
 * @namespace space::randomGen
 * Documentation for space
 */
namespace space::randomGen {

  template<typename Dtype>
  class Uniform {
  private:
    Uniform() : gen_(rd_()), dist_(0, 1) {}

    Uniform(const Uniform &) = default;

    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_real_distribution<Dtype> dist_;
    static std::mutex mutex_;

  public:
    inline static Dtype get(Dtype low = 0, Dtype high = 1) {
      static Uniform singleton;
      return low + (high - low) * singleton.dist_(singleton.gen_);
    }

    inline static Dtype lock_get(Dtype low = 0, Dtype high = 1) {
      static Uniform singleton;

      std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
      lock.lock();
      auto random = singleton.dist_(singleton.gen_);
      lock.unlock();

      return low + (high - low) * random;
    }

    template<typename RandGen>
    inline static Dtype get(RandGen &gen, Dtype low = 0, Dtype high = 1) {
      static Uniform singleton;
      return low + (high - low) * singleton.dist_(gen);
    }
  };

  template<typename Dtype>
  std::mutex Uniform<Dtype>::mutex_;

  template<typename Dtype>
  class Logarithm {
  private:
    Logarithm() : gen_(rd_()), dist_(0, 1) {}

    Logarithm(const Logarithm &) = default;

    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_real_distribution<Dtype> dist_;
    static std::mutex mutex_;

  public:
    inline static Dtype get(Dtype low, Dtype high) {
      static Logarithm singleton;
      Dtype log_low = log10(low);
      Dtype log_high = log10(high);
      return pow(10, log_low + (log_high - log_low) * singleton.dist_(singleton.gen_));
    }

    inline static Dtype lock_get(Dtype low, Dtype high) {
      static Logarithm singleton;
      Dtype log_low = log10(low);
      Dtype log_high = log10(high);

      std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
      lock.lock();
      auto random = singleton.dist_(singleton.gen_);
      lock.unlock();

      return pow(10, log_low + (log_high - log_low) * random);
    }

    template<typename RandGen>
    inline static Dtype get(RandGen &gen, Dtype low, Dtype high) {
      static Logarithm singleton;
      Dtype log_low = log10(low);
      Dtype log_high = log10(high);
      return pow(10, log_low + (log_high - log_low) * singleton.dist_(gen));
    }
  };

  template<typename Dtype>
  std::mutex Logarithm<Dtype>::mutex_;

  template<typename Dtype>
  class PowerLaw {
  private:
    PowerLaw() : gen_(rd_()), dist_(0, 1) {}

    PowerLaw(const PowerLaw &) = default;

    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_real_distribution<Dtype> dist_;
    static std::mutex mutex_;

  public:
    inline static Dtype get(Dtype alpha, Dtype low, Dtype high) {
      static PowerLaw singleton;
      if (!space::iseq(alpha, -1.0)) {
        auto beta = alpha + 1;
        auto f_low = pow(low, beta);
        auto f_high = pow(high, beta);
        return pow(f_low + (f_high - f_low) * singleton.dist_(singleton.gen_), 1.0 / beta);
      } else {
        return Logarithm<Dtype>::get(low, high);
      }
    }

    inline static Dtype lock_get(Dtype alpha, Dtype low, Dtype high) {
      static PowerLaw singleton;
      if (!space::iseq(alpha, -1.0)) {
        auto beta = alpha + 1;
        auto f_low = pow(low, beta);
        auto f_high = pow(high, beta);

        std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
        lock.lock();
        auto random = singleton.dist_(singleton.gen_);
        lock.unlock();

        return pow(f_low + (f_high - f_low) * random, 1.0 / beta);
      } else {
        return Logarithm<Dtype>::lock_get(low, high);
      }
    }

    template<typename RandGen>
    inline static Dtype get(RandGen &gen, Dtype alpha, Dtype low, Dtype high) {
      static PowerLaw singleton;
      if (!space::iseq(alpha, -1.0)) {
        auto beta = alpha + 1;
        auto f_low = pow(low, beta);
        auto f_high = pow(high, beta);
        return pow(f_low + (f_high - f_low) * singleton.dist_(gen), 1.0 / beta);
      } else {
        return Logarithm<Dtype>::get(gen, low, high);
      }
    }
  };

  template<typename Dtype>
  std::mutex PowerLaw<Dtype>::mutex_;

  template<typename Dtype>
  class Normal {
  private:
    Normal() : gen_(rd_()), dist_(0, 1) {}

    Normal(const Normal &) = default;

    std::random_device rd_;
    std::mt19937 gen_;
    std::normal_distribution<Dtype> dist_;
    static std::mutex mutex_;

  public:
    inline static Dtype get(Dtype mean = 0, Dtype sigma = 1) {
      static Normal singleton;
      return mean + sigma * singleton.dist_(singleton.gen_);
    }

    inline static Dtype lock_get(Dtype mean = 0, Dtype sigma = 1) {
      static Normal singleton;

      std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
      lock.lock();
      auto random = singleton.dist_(singleton.gen_);
      lock.unlock();

      return mean + sigma * random;
    }

    template<typename RandGen>
    inline static Dtype get(RandGen &gen, Dtype mean = 0, Dtype sigma = 1) {
      static Normal singleton;
      return mean + sigma * singleton.dist_(gen);
    }
  };

  template<typename Dtype>
  std::mutex Normal<Dtype>::mutex_;

  template<typename Dtype>
  class Maxwell {
  private:
    Maxwell() : gen_(rd_()), dist_(0, 1) {}

    Maxwell(const Maxwell &) = default;

    std::random_device rd_;
    std::mt19937 gen_;
    std::normal_distribution<Dtype> dist_;
    static std::mutex mutex_;

  public:
    inline static Dtype get(Dtype sigma = 1) {
      static Maxwell singleton;
      auto x = singleton.dist_(singleton.gen_);
      auto y = singleton.dist_(singleton.gen_);
      auto z = singleton.dist_(singleton.gen_);

      return sigma * sqrt(x * x + y * y + z * z);
    }

    inline static Dtype lock_get(Dtype sigma = 1) {
      static Maxwell singleton;

      std::unique_lock<std::mutex> lock(mutex_, std::defer_lock);
      lock.lock();
      auto x = singleton.dist_(singleton.gen_);
      auto y = singleton.dist_(singleton.gen_);
      auto z = singleton.dist_(singleton.gen_);
      lock.unlock();

      return sigma * sqrt(x * x + y * y + z * z);
    }

    template<typename RandGen>
    inline static Dtype get(RandGen &gen, Dtype sigma = 1) {
      static Maxwell singleton;
      auto x = singleton.dist_(gen);
      auto y = singleton.dist_(gen);
      auto z = singleton.dist_(gen);

      return sigma * sqrt(x * x + y * y + z * z);
    }
  };

  template<typename Dtype>
  std::mutex Maxwell<Dtype>::mutex_;
}  // namespace space::randomGen

namespace space::random {
  using Uniform = space::randomGen::Uniform<double>;
  using Logarithm = space::randomGen::Logarithm<double>;
  using PowerLaw = space::randomGen::PowerLaw<double>;
  using Normal = space::randomGen::Normal<double>;
  using Maxwell = space::randomGen::Maxwell<double>;
}  // namespace space::random

#endif  //SPACEHUB_RAND_GENERATOR_HPP
