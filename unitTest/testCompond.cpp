#include "../rdFloat.h"
#include <stdio.h>
#include <chrono>
using namespace std::chrono;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> resolutionClock;
int main(int argc, char**argv)
{
    resolutionClock start;
    resolutionClock now;
    
    
    ddouble a = 1, b = 2, c = 2;
    double e = 2, f =0;
    int m =3;
    ddouble sum = 0;
    start         = high_resolution_clock::now();
    for(ddouble i = 1 ; i < 100000000; i+=1)
    {
        sum += 1.0/(i+e);
    }
    now           = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(now - start);
    auto runtime  = double(duration.count()) * milliseconds::period::num / milliseconds::period::den;
    printf("%lf\n",runtime);
    printf("sum=%lf %le\n",sum.num,sum.err);
    if(b==c)
        printf("!!\n");
    if(f==0)
        printf("!!!\n");
    return 0;
    
}
