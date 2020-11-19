//
// Created by yihan on 4/4/19.
//

#pragma once

#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace space::tools {

    static std::string auto_name(std::string const &prefix = "space_") {
        static int duplicate = 1;
        static std::string last_name;

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%d-%m-%Y_%H:%M:%S");

        auto name = prefix + oss.str();

        if (name == last_name) {
            return name + '(' + std::to_string(duplicate++) + ").dat";
        } else {
            duplicate = 1;
            last_name = name;
            return name + ".dat";
        }
    }

    static inline std::string make_name(std::string const &prefix, int suffix, std::string const &extension = ".dat") {
        return prefix + std::to_string(suffix) + extension;
    }
}  // namespace space::tools
