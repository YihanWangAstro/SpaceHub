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
 * @file IO.hpp
 *
 * Header file.
 */
#pragma once

#include <array>
#include <iostream>
#include <tuple>
#include <vector>

namespace space {

    /**
     * @brief Print variables to an output stream
     *
     * @tparam Ostream Output stream.
     * @tparam Args Variadic type(any).
     * @param[in,out] os Output stream.
     * @param[in] args Variables.
     * @return Output stream.
     */
    template <typename Ostream, typename... Args>
    auto &print(Ostream &os, Args &&... args) {
        (os << ... << std::forward<Args>(args));
        return os;
    }

    /**
     * @brief Input variables from an input stream
     *
     * @tparam Istream Input stream.
     * @tparam Args Variadic type(any).
     * @param[in,out] istream Input stream.
     * @param[out] args Variables.
     * @return Input stream.
     */
    template <typename Istream, typename... Args>
    auto &input(Istream &istream, Args &&... args) {
        (istream >> ... >> std::forward<Args>(args));
        return istream;
    }

    /**
     * @brief Print variables to std::cout.
     *
     * @tparam Args Variadic type(any).
     * @param[in] args Variables.
     * @return std::cout.
     */
    template <typename... Args>
    auto &std_print(Args &&... args) {
        (std::cout << ... << std::forward<Args>(args));
        return std::cout;
    }

    /**
     * @brief Input variables from std::cin.
     *
     * @tparam Args Variadic type(any).
     * @param[out] args Variables.
     * @return std::cin.
     */
    template <typename... Args>
    auto &std_input(Args &&... args) {
        (std::cin >> ... >> std::forward<Args>(args));
        return std::cin;
    }

    /**
     * @brief Print comma seperated value to an ostream.
     *
     * @tparam Ostream Type of ostream.
     * @tparam Arg Variadic type(any).
     * @tparam Args Variadic type(any).
     * @param[out] out output stream.
     * @param[in] arg Variable.
     * @param[in] args Variables.
     * @return auto& output stream.
     */
    template <typename Ostream, typename Arg, typename... Args>
    auto &print_csv(Ostream &out, Arg &&arg, Args &&... args) {
        out << arg;
        (..., (out << ',' << std::forward<Args>(args)));
        return out;
    }

    /**
     * @brief Print space seperated value to an ostream.
     *
     * @tparam Ostream Type of ostream.
     * @tparam Args Variadic type(any).
     * @param[out] out Output stream
     * @param[in] args Variables.
     */
    template <typename Ostream, typename... Args>
    void display(Ostream &out, Args &&... args) {
        (..., (out << std::forward<Args>(args) << ' '));
    }

    template <typename... Args>
    std::ostream &operator<<(std::ostream &out, std::tuple<Args...> const &tup) {
        std::apply(
            [&](auto &&arg, auto &&... args) {
                out << arg;
                (..., (out << ' ' << args));
            },
            tup);
        return out;
    }

    template <typename... Args>
    std::istream &operator>>(std::istream &in, std::tuple<Args...> &&tup) {
        std::apply([&](auto &&... args) { (..., (in >> args)); }, tup);
        return in;
    }

    template <typename T, size_t N>
    std::ostream &operator<<(std::ostream &os, std::array<T, N> const &container) {
        for (auto const &c : container) {
            os << c << ' ';
        }
        return os;
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &os, std::vector<T> const &container) {
        for (auto const &c : container) {
            os << c << ' ';
        }
        return os;
    }
}  // namespace space
