
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
                if(file_ == nullptr){
                    printf( "fail to open the file!\n");
                    exit(0);
                }
            }
            ~SynFile(){
                if(file_ == nullptr)
                    fclose(file_);
            }

            template <typename ...Args>
            void write (char const *format, Args&&...args) {
                std::lock_guard<std::mutex> lock(write_mutex_);
                fprintf(file_, format, std::forward<Args>(args)...);
                fflush(file_);
            }

            template <typename ...Args>
            bool read (char const *format, Args&&...args) {
                std::lock_guard<std::mutex> lock(write_mutex_);
                fscanf(file_, format, std::forward<Args>(args)...);
                
                if(feof(file_))
                    return true;
                else
                    return false;
            }

        private:
            FILE* file_{nullptr};
            std::mutex write_mutex_;
        };

        class FileHolder {
        public:
            FileHolder (std::shared_ptr<SynFile> synfile) : synfile_(synfile) {}

            template <typename ...Args>
            void write (char const *format, Args&&...args) {
                synfile_->write(format, std::forward<Args>(args)...);
            }

            template <typename ...Args>
            bool read (char const *format, Args&&...args) {
                return synfile_->read(format, std::forward<Args>(args)...);
            }
        private:
            std::shared_ptr<SynFile> synfile_;
        };
    }
}

#endif
