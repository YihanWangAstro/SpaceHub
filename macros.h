
#ifndef MACROS_H
#define MACROS_H
#include <algorithm> 
#include <memory>

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

constexpr double PI        = 3.14159265358979323;
constexpr double AU        = (PI/648000);
constexpr double PC        = 1;
constexpr double M_SOLAR   = 1;
constexpr double M_JUPITER = 0.9547919E-3;
constexpr double R_SOLAR   = 2.25461E-8;
constexpr double YEAR      = 6.694685210039141E-08;
constexpr double DAY       = YEAR/365.25636042;
constexpr double G         = 1;
constexpr double V_UNIT    = 6.54589713446219E-2;
constexpr double C         = 299792.458/V_UNIT;
constexpr double KM        = 3.2407557442395564e-14;

enum       PARTICTYPE     { NEUTRONSTAR, STAR, BLACKHOLE, POsize_t, NONE = 0 };
enum       EVENTTYPE      { TDE, MERGE, ESCAPE, DISRUPTED, UNEVENTFUL, HVS };
enum class INTEGRATORTYPE { DKDLEAPFROG, KDKLEAPFROG, SYM4, PEFRL, SYM6, SYM8, SYM10 };
enum class SYSTEMTYPE     { PLAIN, CHAIN };
enum class REGUTYPE       { LOGH, TTL, NONE };
enum class ITERTYPE       { BSITER, SEQITER };
#endif
