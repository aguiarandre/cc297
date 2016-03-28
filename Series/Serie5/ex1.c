#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


#define LAMBDA -0.97736

double rk4(double, double);
double rk2(double, double);
double rk1(double, double);

int main(int argc, char* argv[])
{
    
    int N_PASSOS;
    int scheme;
    
    printf("Escolha o método: \n");
    printf("1 para RK-4 \n");
    printf("2 para RK-2 \n");
    printf("3 para RK-1 \n");
    scanf("%d", &scheme);
    
  //  printf("Escolha o número de passos: \n");
  //  scanf("%d", &N_PASSOS);

    double T_MAX = 2.0;
    double w;                    /** Solução */
    double h;                    /** DELTA T */
    double w0 = 1.0;             /** Condição Inicial */
    int n;
    

    switch(scheme)
    {
    case 1:
    
    /** RK-4, item B */
    
    /** 5 Passos no tempo */
    N_PASSOS = 5;
    w = w0;
    n = 0;
    h = T_MAX / N_PASSOS;
    while (n < N_PASSOS)
    {
        w = rk4( h, w );
        n++;
    }
    printf("RK4, N_PASSOS = %02.d,  w = %.11f\n", N_PASSOS, w);

    /** 3 Passos no tempo */
    N_PASSOS = 3;
    w = w0;
    n = 0;
    h = T_MAX / N_PASSOS;
    
    while (n < N_PASSOS)
    {
        w = rk4( h, w );
        n++;
    }
    printf("RK4, N_PASSOS = %02.d,  w = %.11f\n", N_PASSOS, w);

    /** 10 Passos no tempo */

    N_PASSOS = 10;
    w = w0;
    n = 0;
    h = T_MAX / N_PASSOS;

    while (n < N_PASSOS)
    {
        w = rk4( h, w );
        n++;
    }
    printf("RK4, N_PASSOS = %02.d,  w = %.11f\n", N_PASSOS, w);
    break;
    
    case 2:
    
    N_PASSOS = 10;
    w = w0;
    n = 0;
    h = T_MAX / N_PASSOS;
    
    while (n < N_PASSOS)
    {
        w = rk2( h, w );
        n++;
    }
    printf("RK2, N_PASSOS = %02.d,  w = %.11f\n", N_PASSOS, w);

    break;
    
    case 3:
        
    N_PASSOS = 20;
    w = w0;
    n = 0;
    h = T_MAX / N_PASSOS;
    
    while (n < N_PASSOS)
    {
        w = rk1( h, w );
        n++;
    }
    printf("RK1, N_PASSOS = %02.d,  w = %.11f\n", N_PASSOS, w);

    break;
    
    default:
    
    printf("Opção errada!\n");
    break;
        
    }
    
    return 0;
}

double rk4(double h, double w)
{
    
    double wTil, wBar, wChapeu;
    
    wTil = w + (0.5 * h * LAMBDA * w);
    wBar = w + (0.5 * h * LAMBDA * wTil);
    wChapeu = w + h * LAMBDA * wBar;
    
    w = w + (1.0/6.0) * h * LAMBDA * (wChapeu + 2 * (wBar + wTil) + w);
    return w;
}


double rk2(double h, double w)
{
    
    double wTil;
    
    wTil = w + (0.5 * h * LAMBDA * w);

    w = w +  h * LAMBDA *  wTil;
    
    return w;
}

double rk1(double h, double w)
{
    
    w = w +  h * LAMBDA *  w;
    
    return w;
}