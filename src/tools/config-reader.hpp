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
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file config-reader.hpp
 *
 * Header file.
 */
#pragma once

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "../dev-tools.hpp"

/**
 * @namespace space::tools
 * space name for tools
 */
namespace space::tools {

    /**
     *
     */
    enum class ConfigDtype { Integer, Float, String, Empty };
    /*---------------------------------------------------------------------------*\
          Class ConfigReader Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     */
    class ConfigReader {
       public:
        explicit ConfigReader(std::string const &file_name, char divider = '=', char commenter = '#');

        explicit ConfigReader(char const *file_name, char divider = '=', char commenter = '#');
        ;

        template <typename T>
        T get(std::string const &key);

       private:
        std::unordered_map<std::string, std::string> map_;
    };

    /*---------------------------------------------------------------------------*\
          Class ConfigReader Implementation
    \*---------------------------------------------------------------------------*/
    ConfigReader::ConfigReader(std::string const &file_name, char divider, char commenter) {
        std::fstream file(file_name);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

                if (line[0] == commenter || line.empty()) continue;

                auto split = line.find(divider);
                auto key = line.substr(0, split);
                auto val = line.substr(split + 1);
                map_[key] = val;
            }
        } else {
            spacehub_abort("Cannot open the configure file: ", file_name, "\r\n");
        }
    }

    ConfigReader::ConfigReader(char const *file_name, char divider, char commenter)
        : ConfigReader(std::string(file_name), divider, commenter) {}

    template <typename T>
    T ConfigReader::get(std::string const &key) {
        auto iskey = (map_.end() != map_.find(key));
        if (iskey) {
            std::stringstream ss(map_[key]);
            T value;
            ss >> value;
            return value;
        } else {
            spacehub_abort("Invalid key for configure file!");
        }
    }

    template <typename... Args>
    void read_command_line(int argc, char **argv, Args &&...args) {
        constexpr size_t n = sizeof...(Args);

        if (argc != n + 1) {
            spacehub_abort("Wrong args number!");
        } else {
            std::stringstream ss;
            for (int i = 1; i < argc; ++i) {
                ss << argv[i] << ' ';
            }
            ((ss >> std::forward<Args>(args)), ...);
        }
    }

    class OptArg {
       public:
        inline operator float() { return std::stof(opt_); }

        inline operator double() { return std::stod(opt_); }

        inline operator long double() { return std::stold(opt_); }

        inline operator int() { return std::stoi(opt_); }

        inline operator long() { return std::stol(opt_); }

        inline operator long long() { return std::stoll(opt_); }

        inline operator unsigned long() { return std::stoul(opt_); }

        inline operator unsigned long long() { return std::stoull(opt_); }

        inline operator char const *() { return opt_.c_str(); }

        inline operator std::string() { return opt_; }

        inline operator bool() { return static_cast<bool>(std::stoi(opt_)); }

        template <typename T>
        T as() {
            if constexpr (std::is_same_v<T, float>) {
                return std::stof(opt_);
            } else if constexpr (std::is_same_v<T, double>) {
                return std::stod(opt_);
            } else if constexpr (std::is_same_v<T, long double>) {
                return std::stold(opt_);
            } else if constexpr (std::is_same_v<T, int>) {
                return std::stoi(opt_);
            } else if constexpr (std::is_same_v<T, long>) {
                return std::stol(opt_);
            } else if constexpr (std::is_same_v<T, long long>) {
                return std::stoll(opt_);
            } else if constexpr (std::is_same_v<T, unsigned long>) {
                return std::stoul(opt_);
            } else if constexpr (std::is_same_v<T, unsigned long long>) {
                return std::stoull(opt_);
            } else if constexpr (std::is_same_v<T, char *>) {
                return opt_.c_str();
            } else if constexpr (std::is_same_v<T, std::string>) {
                return opt_;
            }
        }

        friend std::istream &operator>>(std::istream &is, OptArg &arg) {
            is >> arg.opt_;
            return is;
        }

        friend std::ostream &operator<<(std::ostream &os, OptArg const &arg) {
            os << arg.opt_;
            return os;
        }

       private:
        std::string opt_;
    };

#define READ_CMD_LINE(ARGC, ARGV, ...) \
    tools::OptArg __VA_ARGS__;         \
    read_command_line(ARGC, ARGV, __VA_ARGS__);

}  // namespace space::tools
