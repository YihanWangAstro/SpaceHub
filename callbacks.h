
#ifndef CALLBACKS_H
#define CALLBACKS_H

#include "devTools.h"
#include "tools/timer.h"
#include <fstream>
#include <iomanip>

namespace SpaceH {

    namespace CallBack {

        template<typename ParticleSys>
        class DefaultWriter {
        public:
            using Scalar = typename ParticleSys::Scalar;

            DefaultWriter(const char *path, Scalar end_time, size_t output_num = 5000)
                    : write_interval_(end_time / output_num){

                os = std::make_shared<std::ofstream>(path);
                if (os->is_open()) {
                    (*os) << std::scientific << std::setprecision(16);
                } else {
                    ERR_MSG(("Fail to open the file" + std::string(path)).c_str());
                }
            }

            inline void operator()(ParticleSys &partc) {
                if (partc.time() >= write_time_) {
                    (*os) << partc;
                    write_time_ += write_interval_;
                }
            }

        private:
            std::shared_ptr<std::ofstream> os;
            Scalar write_time_{0};
            Scalar write_interval_{0};
            size_t step_{0};
        };

        template<typename ParticleSys>
        class ShowProgressBar{
        public:
            using Scalar = typename ParticleSys::Scalar;
            explicit ShowProgressBar(Scalar end_time) : bar(end_time) {}

            inline void operator()(ParticleSys &partc) {
                bar.autoShow(partc.time());
            }
        private:
          Progress bar;
        };

    }
}
#endif

