//
// Created by yihan on 2/15/19.
//

#ifndef SPACEHUB_CONFIG_READER_H
#define SPACEHUB_CONFIG_READER_H

#include <fstream>
#include <unordered_map>
#include <algorithm>
#include "dev_tools.h"

namespace SpaceH {
    class ConfigReader : public std::unordered_map<std::string, double> {
    public:
        ConfigReader(std::string const &file_name, char divider = '=', char commenter = '#') {
            std::fstream file(file_name);
            if (file.is_open()) {
                std::string line;
                while (std::getline(file, line)) {
                    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

                    if (line[0] == commenter || line.empty())
                        continue;

                    auto devider_pos = line.find(divider);
                    auto key = line.substr(0, devider_pos);
                    auto value = atof(line.substr(devider_pos + 1).c_str());
                    (*this)[key] = value;
                }
            } else {
                SPACEHUB_ABORT("Cannot open the configure file: ", file_name, "\r\n");
            }
        }

        explicit ConfigReader(char const *file_name, char divider = '=', char commenter = '#')
                : ConfigReader(std::string(file_name), divider, commenter) {};

        friend std::ostream &operator<<(std::ostream &os, ConfigReader const &config) {
            for (auto&[key, value] : config) {
                SpaceH::print(os, key, '=', value, '\n');
            }
        }
    };
}


#endif //SPACEHUB_CONFIG_READER_H
