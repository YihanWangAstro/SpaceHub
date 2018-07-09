#ifndef TIMMER_H
#define TIMMER_H
#include <chrono>

class ProgressBar
{
public:
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::time_point<Clock> ClockTime;
    typedef std::chrono::milliseconds ms;

    ProgressBar() = delete;
    ProgressBar(double upperLimit, unsigned int precision = 10000) : end(upperLimit), step(upperLimit/precision), current(0.0)
    {
        start();
    }
    
    void autoShow(double time)
    {
        if(time >= current)
        {
            auto now       = Clock::now();
            auto duration  = duration_cast<ms>(now - startTime);
            double percent = time/end;
            runtime        = double(duration.count()) * ms::period::num / ms::period::den;
            
            NOTICE::EraseLine();
            printf("Computing >>>" ANSI_COLOR_RED " %5.1lf%% " ANSI_COLOR_RESET "     Time rest >>>"
                   ANSI_COLOR_RED " %9.3lf" ANSI_COLOR_RESET "s     Runtime >>>" ANSI_COLOR_RED
                   " %9.3lf" ANSI_COLOR_RESET "s",percentage*100, runtime*(1 -  percentage)/percentage, runtime);
            fflush(stdout);
            current += step;
        }
    }
    
    void reset(double upperLimit, int precision = 10000)
    {
        end       = upperLimit;
        step      = upperLimit/precision;
        current   = 0.0;
        startTime = start();
    }
        
    void start()
    {
        startTime = Clock::now();
    }
    
    double getTime()
    {
        auto now      = Clock::now();
        auto duration = duration_cast<ms>(now - startTime);
        return  double(duration.count()) * ms::period::num / ms::period::den;
    }
    
private:
    ClockTime startTime;
    double    end;
    double    step;
    double    current;
};
#endif
