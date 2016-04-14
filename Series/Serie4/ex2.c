/**
 * Exercicio 2 - Serie 4
 *
 */ 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846
#define T_MAX 2

void multiply(double* a, double* b, double* c, int );
double method(double left, double mid, double right, double dx);

int main(int argc, char* argv[])
{
    if (argc > 2)
    {
        printf("Usage: ./ex2 \n");
    }
    
    /** Declarações ! */
    
    int M = 5;    
    int n = M+1;
    int N_PASSOS = 5;
    double h = (double) T_MAX/N_PASSOS; /** Time-Step */
    int nPassos = 0;
    double sum = 0;
    
    double * x = malloc ( n * sizeof(double));
    double * lambda = malloc( n * sizeof(double));
    double * u0_real = malloc( (n)*sizeof(double) );
    double * u = malloc( (n)*sizeof(double) );    /** u - physical space */    
    double * w = malloc( (n)*sizeof(double) );
    double * uBar = malloc( (n)*sizeof(double) );
    double * uTil = malloc( (n)*sizeof(double) );
    double * uHat = malloc( (n)*sizeof(double) );
    double * rhs = malloc( (n)*sizeof(double) );
    double * rhsTil = malloc( (n)*sizeof(double) );
    double * rhsHat = malloc( (n)*sizeof(double) );
    double * rhsBar = malloc( (n)*sizeof(double) );
    double ** autovec = malloc( n * sizeof(double*));
    double ** invAutovec = malloc( n * sizeof(double*));
    for (int i = 0; i < n ; i++)
    {
        autovec[i] = malloc( n * sizeof(double));
        invAutovec[i] = malloc( n * sizeof(double));
    }
   
    double w0_onda[] = { 1.0, 0.1, 0.1, 0.1, 0.1 }; /** Initial Condition on wave-space */
    x[0] = 0.0;                                     /** Left boundary of the mesh*/
    double dx = PI / ( M + 1);                      /** Mesh spacing */
    
   /** Mesh generation ! */
   
    for (int i = 0, m = M+2; i < m; i++)
    {
        x[i] = x[0] + (dx * i); 
    }

    /** Eigen-vectors/values (and its inverse) calculations */
    
    for ( int m = 1; m < n; m++)
    {
        lambda[m-1] = - (2.0/(dx*dx)) * (1.0 - cos( ( (double) (m * PI)/(M + 1) )) );
        
        for(int j = 1; j < n ; j++)
        {
            autovec[m-1][j-1] = sin( (j * m * PI)/ (M+1) ) ;
            invAutovec[m-1][j-1] = (2.0/(M+1)) * autovec[m-1][j-1];
        }
    }
    
 
    /** Cálculo de u0 no espaço físico */

    for (int j = 1 ; j < n ; j++)
    {
        for ( int m = 1 ; m < n ; m++)
        {
            sum = sum + autovec[m-1][j-1] * w0_onda[m-1];
        }
        u[j] = sum;
        u0_real[j] = sum;
        sum = 0;
    }
    

    /** Condição de contorno! */
    u[0] = 0.0;
    u[n] = 0.0;
    u0_real[0] = 0.0;
    u0_real[n] = 0.0;

    /** 
     * LETRA A
     */
     
    while ( nPassos < N_PASSOS)
    {
    
        for (int i = 1; i < n; i++) /** Só há atualização dos pontos internos! */
        {
            rhs[i] = method(u[i-1],u[i],u[i+1],dx);
            uTil[i] = u[i] + ( 0.5 * h * rhs[i] );
            
        }
        for (int i = 1; i < n; i++)
        {
            rhsTil[i] = method(uTil[i-1], uTil[i], uTil[i+1],dx);
            uBar[i] = u[i] + (0.5 * h * rhsTil[i]);
        }
        
        for (int i = 1; i < n; i++) 
        {
            rhsBar[i] = method(uBar[i-1], uBar[i], uBar[i+1],dx);
            uHat[i] = u[i] + ( h * rhsBar[i]);
        }
        
        for (int i = 1; i < n; i++) 
        {
            rhsHat[i] = method(uHat[i-1], uHat[i], uHat[i+1],dx);
      
            u[i] = u[i] + (1.0/6.0) * h * ( rhsHat[i] + 2*(rhsBar[i] + rhsTil[i]) + rhs[i]);
        }
    
        nPassos++;        
        
    }
    
    sum = 0;
    double w0[M+1];
    
    for (int j = 1 ; j < n ; j++)
    {
        for ( int m = 1 ; m < n ; m++)
        {
            sum = sum + invAutovec[m-1][j-1] * u[m];
            w0[m] = w0_onda[m-1];
             
        }
    
        w[j] = sum;
        
        sum = 0;
    
    }
    
    w[0] = 0.0;
    w[n] = 0.0;
    w0[0] = 0.0;
    w0[n] = 0.0;
     
     
        printf("      x         u0             u          w0            w\n");
        for (int i = 0; i < n+1; i++) 
        {
            
            printf("%0.10f, %0.10f, %0.10f %0.10f %0.10f\n", x[i], u0_real[i], u[i], w0[i], w[i] );
        } 
        printf("\n");
       
        for (int i = 1; i < n; i++) /** Só há atualização dos pontos internos! */
        {
            
            printf("lambda[%d] = %f     ", i-1, lambda[i-1] );
        } 
     printf("\n h = %f, maxH = %f\n", h, -2.78529356/lambda[M-1]);   
    
    N_PASSOS = 10;
    h = (double) T_MAX/N_PASSOS; /** Time-Step */
    
    /** Item D - Runge Kutta 2 */
    printf("RK-2 \n");
        
        
    /** Cálculo de u0 no espaço físico */

    for (int j = 1 ; j < n ; j++)
    {
        for ( int m = 1 ; m < n ; m++)
        {
            sum = sum + autovec[m-1][j-1] * w0_onda[m-1];
        }
        u[j] = sum;
        u0_real[j] = sum;
        sum = 0;
    }
     
     
    nPassos = 0;       
        
    while ( nPassos < N_PASSOS)
    {
    
        for (int i = 1; i < n; i++) /** Só há atualização dos pontos internos! */
        {
            rhs[i] = method(u[i-1],u[i],u[i+1],dx);
            uTil[i] = u[i] + ( 0.5 * h * rhs[i] );
            
        }
        for (int i = 1; i < n; i++)
        {
            rhsTil[i] = method(uTil[i-1], uTil[i], uTil[i+1],dx);
            u[i] = u[i] + ( h * rhsTil[i]);
        }
        
        nPassos++;        
        
    }  
    
    sum = 0;
    
    for (int j = 1 ; j < n ; j++)
    {
        for ( int m = 1 ; m < n ; m++)
        {
            sum = sum + invAutovec[m-1][j-1] * u[m];

        }
    
        w[j] = sum;
        
        sum = 0;
    
    }
    

        printf("      x         u0             u          w0            w\n");
        for (int i = 0; i < n+1; i++) 
        {
            
            printf("%0.10f, %0.10f, %0.10f %0.10f %0.10f\n", x[i], u0_real[i], u[i], w0[i], w[i] );
        } 
      printf("\n");
       
        for (int i = 1; i < n; i++) /** Só há atualização dos pontos internos! */
        {
            
            printf("lambda[%d] = %f     ", i-1, lambda[i-1] );
        }    
        
        
        
  printf("\n h = %f, maxH = %f\n", h, -2.0/lambda[M-1]);  
        

    for (int i = 0; i < n; i++)
    {
        free(autovec[i]);
        free(invAutovec[i]);
    }
    
    free(autovec);
    free(invAutovec);
    free(x);
    free(u);
    free(w);
    free(u0_real);
    free(uHat);
    free(uBar);
    free(uTil);
    free(rhs);
    free(rhsHat);
    free(rhsBar);
    free(rhsTil);
    free(lambda);
    
    
    
    return 0;
}


  /*  
    double wTil, wBar, wChapeu;
    wTil = w + (0.5 * h * RHS);
    wBar = w + (0.5 * h * LAMBDA * wTil);
    wChapeu = w + h * LAMBDA * wBar;
    
    w = w + (1.0/6.0) * h * LAMBDA * (wChapeu + 2 * (wBar + wTil) + w); 
    */


  /*
    double aux;
    
    for (int k = 0; k < M; k++ )
    {   
        for( int i = 0; i < M ; i ++)
        {
            for ( int j = 0; j < M ; j++)
            {
                aux = aux + autovec[k][j]*invAutovec[j][i];
            }
        
            identity[k][i] = aux;
           // printf("identity[%d][%d] = %f\n", k, i, identity[k][i]);

            aux = 0.0;
            
        }
        
    }
    */

double method(double left, double mid, double right, double dx)
{
    double ni = 1.0;
    return (ni/(dx*dx) ) * (left - 2*mid + right);
}

