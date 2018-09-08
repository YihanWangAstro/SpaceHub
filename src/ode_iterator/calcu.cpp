#include<iostream>
#include<math.h>
long double velTaylorCoef(long double h, size_t index)//make it const later
{
    double var = 1;
    size_t order = index + 2;
    /*for(size_t i = 0 ; i < order; ++i)
        var*=h;
    
    return var/order;*/
    return pow(h,order)/order;
}

long double posTaylorCoef(long double h, size_t index)//make it const later
{
    double var = 1;
    size_t order = index + 3;
    /*for(size_t i = 0 ; i < order; ++i)
        var*=h;
    
    return var/(order*(order-1));*/
    return pow(h,order)/(order*(order-1));
}

int main()
{
    long double H[9] = {0,0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1};
    
    
    
    /*for(size_t i = 0 ; i < 8; ++i)
    {
        //printf("%.17Le, ",H[i]);
        for(int j = -1 ; j < 7; ++j)
        {
            printf("%.17Le, ", velTaylorCoef(H[i],j));
        }
        printf("\n");
    }*/
    //std::cout << std::scientific << std::setprecision(17);
    /*double c[9][9];
    double d[9][9];
    
    for(int i = 0 ; i < 8 ; ++i)
    {
        for(int j = 0 ;j<8 ; ++j)
        {
            printf("%d %d\n",i,j);
            if(i==j)
            {
                c[i][j]=1;
            }
            
            else if(j>i)
                c[i][j]= 0;
            else if(j == 0)
                c[i][j]= -H[i]*c[i-1][0];
            
            else
                c[i][j]=c[i-1][j-1]-H[i]*c[i-1][j];
        }
    }
    
    for(int i = 0 ; i < 9; ++i)
    {
        for(int j = 0 ; j < 9; ++j)
        {
            printf("%.17le, ",c[j][i]);
        }
        printf("\n");
    }
    
    
    for(int i = 0 ; i < 8 ; ++i)
    {
        for(int j = 0 ;j<8 ; ++j)
        {
            printf("%d %d\n",i,j);
            if(i==j)
            {
                d[i][j]=1;
            }
            else if(j>i)
                d[i][j]= 0;
            else if(j == 0)
                d[i][j]= H[1]*d[i-1][0];
            
            else
                d[i][j]=d[i-1][j-1] + H[j+1]*d[i-1][j];
        }
    }
    
    for(int i = 0 ; i < 9; ++i)
    {
        for(int j = 0 ; j < 9; ++j)
        {
            printf("%.17le, ",d[j][i]);
        }
        printf("\n");
    }*/
    
    long double r[9][9];
    long double rr[9][9];
    for(int i = 0 ; i < 8 ; ++i)
    {
        for(int j = 0 ; j<8 ; ++j)
        {
            printf("%d %d\n",i,j);
            if(j == 0)
            {
                r[i][j]=1/H[i];
                rr[i][j] = H[i];
            }
            else if(j>i+1)
            {
                r[i][j]=0;
                rr[i][j]=0;
            }
            else
            {
                r[i][j]=1/(H[i]-H[j]);
                rr[i][j]=(H[i]-H[j]);
            }
        }
    }
    
    for(int i = 0 ; i < 9; ++i)
    {
        for(int j = 0 ; j < 9; ++j)
        {
            printf("%.17Le, ",r[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0 ; i < 9; ++i)
    {
        for(int j = 0 ; j < 9; ++j)
        {
            printf("%.17Le, ",rr[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    long double pod = 0;
    for(int i = 0 ; i < 8; ++i)
    {
        for(int j = 0 ; j < 8; ++j)
        {
            pod = 1;
            for(int k = i-1 ; k >= j; k--)
            {
                pod *= r[i][k];
            }
            printf("%.17Le, ",pod);
        }
        printf("\n");
    }
    
    
    
}
