
#ifndef SPACEHUB_TIMER_H
#define SPACEHUB_TIMER_H
#include <chrono>
#include <stdio.h>
#include <unistd.h>


namespace SpaceH
{
    class Timer
    {
        using Time = std::chrono::time_point<std::chrono::high_resolution_clock>;
    public:
        void start()
        {
            active_ = true;
            start_  = std::chrono::high_resolution_clock::now();
        }
        
        float getTime()
        {
            if(active_)
            {
                auto now = std::chrono::high_resolution_clock::now();
                auto len = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_);
                return  duration_ + static_cast<float>(len.count()) * std::chrono::milliseconds::period::num
                                                                    / std::chrono::milliseconds::period::den;
            }
            else
                return duration_;
        }
        
        void pause()
        {
            duration_ = getTime();
            active_   = false;
        }
        
        void reset()
        {
            active_   = false;
            duration_ = 0;
        }
    private:
        Time start_;
        float duration_{0};
        bool active_{false};
    };
    
    class Progress
    {
    public:
        Progress() = delete;
        Progress(float end) : pro_end_(end), den_(500), pro_ratio_(0)
        {
            ch_[0] = '|', ch_[1] = '/',ch_[2] = '-',ch_[3] = '\\';
            for(int i = 0 ; i < 52; i++)
                done_[i] = '\0';
            
            for(int i = 0 ; i < 51; i++)
                left_[i] = ' ';
            left_[50] = '\0';
        }
        
        void start()
        {
            clock.start();
        }
        
        float getTime()
        {
            return clock.getTime();
        }
        
        void autoShow(float current)
        {
            float r = current / pro_end_ ;
            
            if(r * den_ >= pro_ratio_ && pro_ratio_ <= den_)
            {
                float time = clock.getTime();
                int bar_loc_ = r * 50;
                
                printf("\033[0m\033[?25l\033[42m\033[30m");
                printf("[%c]%s%5.1f%% \033[0m %s \033[0m \033[32mRun time : %.3f s \r\033[0m", ch_[bar_loc_%4], done_, r*100, left_,time);
                fflush(stdout);
                done_[bar_loc_] = ' ';
                left_[50-bar_loc_] = '\0';
                pro_ratio_++;
                
                if(current >= pro_end_)
                    printf("\033[32m\033[90C\033[0m\n");
                    
                printf("\033[?25h");
            }
        }
    private:
        char done_[52];
        char left_[51];
        Timer clock;
        char ch_[4];
        float pro_end_;
        int den_;
        int pro_ratio_;
    };
    
}
#endif
