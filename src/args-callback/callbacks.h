
#ifndef CALLBACKS_H
#define CALLBACKS_H


#include "../dev-tools.h"
#include "../tools/timer.h"
#include <fstream>
#include <iomanip>
#include <memory>

namespace space::ArgsCallBack {


    template <typename Ostream>
    class BaseWriter {
    public:
        BaseWriter(Ostream& os_, double start_, double end_, size_t snapshot_num = 5000, bool flush = false)
                : os{os_}, write_time_{start_}, end_time_{end_}, write_interval_{(end_ - start_) / snapshot_num}, flush_immediate_{flush} {}

        template <typename ParticleSys>
        inline void operator()(ParticleSys &ptc) {
            auto t = ptc.time();
            if (t >= write_time_ && write_time_ <= end_time_) {
                os << ptc << '\n';
                if (flush_immediate_) {
                    os << std::endl;
                }
                write_time_ += write_interval_;
            }
        }

        BaseWriter(BaseWriter const &other) :
                os{other.os},
                write_time_{other.write_time_},
                end_time_{other.end_time_},
                write_interval_{other.write_interval_},
                flush_immediate_{other.flush_immediate_} {}

        void reset_output_params(double start_, double end_, size_t snapshot_num = 5000, bool flush = false) {
            write_time_ = start_;
            end_time_ = end_;
            write_interval_ = (end_ - start_) / snapshot_num;
            flush_immediate_ = flush;
        }

        template<typename T>
        friend BaseWriter &operator<<(BaseWriter &wtr, T const &d) {
            wtr.os << d;
            return wtr;
        }

    private:
        Ostream &os;
        double write_time_{0};
        double end_time_{0};
        double write_interval_{0};
        bool flush_immediate_{false};
    };

    class DefaultWriter {
    public:
        DefaultWriter(std::string const& file_name, double start_, double end_, size_t snapshot_num = 5000, bool flush = false) :
                fstream_{std::make_shared<std::ofstream>(file_name)},
                writer_{*fstream_, start_, end_, snapshot_num, flush} {
            if (!fstream_->is_open()) {
                SPACEHUB_ABORT("Fail to open the file " + file_name);
            } else {
                (*fstream_) /*<< std::fixed*/ << std::setprecision(16);
            }
        }
        DefaultWriter(DefaultWriter const&) = default;

        DefaultWriter&operator=(DefaultWriter const&) = default;

        DefaultWriter(DefaultWriter && other) noexcept : fstream_(std::move(other.fstream_)), writer_(std::move(other.writer_)) {}

        template <typename ParticleSys>
        inline void operator()(ParticleSys &ptc) {
            writer_(ptc);
        }

        void reset_output_params(double start_, double end_, size_t snapshot_num = 5000, bool flush = false){
            writer_.reset_output_params(start_, end_, snapshot_num, flush);
        }

        template <typename T>
        friend DefaultWriter&operator<<(DefaultWriter& wtr, T const& d){
            wtr.writer_ << d;
            return wtr;
        }
    private:
        std::shared_ptr<std::ofstream> fstream_;
        BaseWriter<std::ofstream> writer_;
    };
}
#endif

