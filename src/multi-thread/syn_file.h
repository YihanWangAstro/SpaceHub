
#ifndef SPACEHUB_SYN_FILE_H
#define SPACEHUB_SYN_FILE_H

#include <iostream>
#include <string>
#include <stdio.h>
namespace SpaceH{
    namespace MultiThread{

        class SynFile {
        public:
            SynFile (std::string file_name, std::string mode="w") {
                file_ = fopen(file_name.c_str(), mode.c_str());
                if(file_ == nullptr)
                    throw "fail to open the file!\n";
            }
            ~SynFile(){
                if(file_ == nullptr)
                    fclose(file_);
            }

            template <typename ...Args>
            void write (char const *format, Args...args) {
                std::lock_guard<std::mutex> lock(write_mutex_);
                fprintf(file_, format, args...);
            }

            template <typename ...Args>
            void read (char const *format, Args...args) {
                std::lock_guard<std::mutex> lock(write_mutex_);
                fscanf(file_, format, args...);
            }

        private:
            FILE* file_{nullptr};
            std::mutex write_mutex_;
        };

        class FileHolder {
        public:
            FileHolder (std::shared_ptr<SynFile> synfile) : synfile_(std::move(synfile)) {}

            template <typename ...Args>
            void write (char const *format, Args...args) {
                synfile_->write(format, std::forward<Args>(args)...);
            }

            template <typename ...Args>
            void read (char const *format, Args...args) {
                synfile_->read(format, std::forward<Args>(args)...);
            }
        private:
            std::shared_ptr<SynFile> synfile_;
        };
    }
}

#endif
