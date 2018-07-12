#include <chrono>
#include<iostream>
#include<array>
#include "../kahanNumber.h"
using namespace std::chrono;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> resolutionClock;

template<typename T>
struct dtype
{
    template<typename U>
    static typename U::value_type check(typename U::value_type);
    
    template<typename U>
    static U check(U);
    
    using type = decltype(check<T>(0));
    
};

int main()
{
    resolutionClock start;
    resolutionClock now;
    start         = high_resolution_clock::now();
    
    using kh = SpaceH::kahan<double>;
    
    kh a = 0;
    //using kh = double;
    /*
    size_t N = 10000000;
    kh sum = 1;
    for(kh i = 1; i < N ; i = i + 1)
    {
        sum = sum + 1.0/i;
    }
    
    now           = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(now - start);
    auto runtime  = double(duration.count()) * milliseconds::period::num / milliseconds::period::den;
    
    std::cout << sum << "\n";
    std::cout << runtime << std::endl;*/
    
    
    return 0;
}
