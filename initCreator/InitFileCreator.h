#ifndef INITFILECREATOR_H
#define INITFILECREATOR_H

#include <iomanip>
#include <fstream>
#include "devTools.h"
#include "protoType.h"

namespace SpaceH {
    template<typename Particle>
    class InitCreator {
    public:
        /* Typedef */
        using type   = typename Particle::type;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;

        template<typename T, size_t S>
        using Container = typename type::template Container<T, S>;
        /* Typedef */

        /**
         *
         * @param filePath
         * @param time
         */
        void writeToFile(const char *filePath, Scalar time = 0) {
            std::ofstream outFile(filePath);
            if (outFile.is_open()) {
                outFile << std::scientific << std::setprecision(16);
                outFile << '#';
                size_t size = particleNumber();
                outFile << ' ' << size << ' ' << time << "\r\n";

                for (size_t i = 0; i < size; ++i) {
                    outFile << i << ' ' << particles_[i] << "\r\n";
                }
            } else
                SpaceH::ERR_MSG("Fail to open the initial file to write the Particle Data!", __FILE__, __LINE__);

            outFile.close();
        }

        /**
         *
         * @param newParticle
         */
        void addParticle(const Particle &newParticle) {
            particles_.emplace_back(newParticle);
        }

        /**
         *
         * @return
         */
        size_t particleNumber() {
            return particles_.size();
        }

        /**
         *
         * @param i
         * @return
         */
        inline const Particle &operator[](size_t i) const {
            return particles_.at(i);
        }

        /**
         *
         */
        void clear() {
            particles_.clear();
        }

    private:
        Container<Particle, SpaceH::DYNAMICAL> particles_;
    };
}
#endif
