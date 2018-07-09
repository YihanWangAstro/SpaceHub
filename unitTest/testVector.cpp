#include "../genVector3d.h"
#include "../rdFloat.h"
#include <stdio.h>
int main(int argc, char**argv)
{
    Vector3d<double> v1(0,0,1), v2(1,2,3),v3;
    
    v3 = v1 - v2;
   // printf("%lf %lf %lf %lf %lf %lf %lf",v1.x.num,v1.y.num,v1.z.num,v1.norm(),v3.x.num,v3.y.num,v3.z.num);
    printf("%lf %lf %lf %lf %lf %lf %lf",v1.x,v1.y,v1.z,v1.norm(),v3.x,v3.y,v3.z);
    
    return 0;
    
}
