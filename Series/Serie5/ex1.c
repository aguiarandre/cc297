/**
 * Exercicio 1 - Serie 5 
 * 
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define T_MAX 2

void multiply(double* a, double* b, double* c, int );
double method(double left, double mid, double right, double dx);
double rhs(double CFL, double left, double mid, double right);
bool solveTridiag(double* a, double* b, double* c, double* d, double* correct, int N);


int main(int argc, char* argv[])
{
    if (argc > 2)
    {
        printf("Usage: ./ex1 \n");
    }
    
    /** Declarações ! */
    
    int M = 11;    
    int n = M+1;
    int N_PASSOS = 2;
    double h = (double) T_MAX/N_PASSOS; /** Time-Step */
    int nPassos = 0;
    double sum = 0;
    
    double * x = malloc ( n * sizeof(double));
    double * lambda = malloc( n * sizeof(double));
    double * u0_real = malloc( (n)*sizeof(double) );
    double * u = malloc( (n)*sizeof(double) );    /** u - physical space */    
    double * w = malloc( (n)*sizeof(double) );

    double * a = malloc( (n)*sizeof(double) );      
    double * b = malloc( (n)*sizeof(double) );    
    double * c = malloc( (n)*sizeof(double) );    
    double * d = malloc( (n)*sizeof(double) );

    
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
    double w0[M+1];
    double CFL = 0.5 * h / ( dx*dx ) ;

    
   /** Mesh generation ! */
   
    for (int i = 0, m = M+2; i < m; i++)
    {
        x[i] = x[0] + (dx * i); 
    }

    /** Eigen-vectors/values (and its inverse) calculations */
    
    for ( int m = 1; m < n; m++)
    {
        lambda[m-1] = -( 2.0/ (dx*dx)) * (1.0 - cos( ( (double) (m * PI)/(M + 1) )) );
        
        for(int j = 1; j < n ; j++)
        {
            autovec[m-1][j-1] = sin( (j * m * PI)/ (M+1) );
            invAutovec[m-1][j-1] = (2.0/(M+1)) * autovec[m-1][j-1];
        }
    }
    
 
    /** Ex.3 - Runge Kutta 2 */
    printf("RK-2 \n");
        
        
    /** Cálculo de u0 no espaço físico */

    
    for (int j = 1 ; j < n ; j++)
    {
        u[j] = sin(  (double)j * dx );
        u0_real[j] = u[j];
    }
     
     
    nPassos = 0;       
   /** Condição de contorno! */
    u[0] = 0.0;
    u[n] = 0.0;
    u0_real[0] = 0.0;
    u0_real[n] = 0.0;

    
        
    while ( nPassos < N_PASSOS)
    {
    
        for ( int i = 1; i < n; i++)
        {

            /** Generate the matrix elements */
            a[i] = CFL;
            b[i] = - (2.0*CFL + 1);
            c[i] = CFL;
            
            /** Calculate RHS vector */
            d[i] = rhs(CFL, u[i-1], u[i], u[i+1]);
            
        }
        
        /** solve tridiagonal */   
        solveTridiag(a, b, c, d, u, n);
        
        
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
            
          //  printf("lambda[%d] = %f     ", i-1, lambda[i-1] );
        }    
        
        
        //printf("\n h = %f, maxH = %f\n", h, -2.0/lambda[M-5]);
        
        

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
    free(d);
    free(lambda);
    
    
    
    return 0;
}



double method(double left, double mid, double right, double dx)
{

    double ni = 1.0;
    return (ni/(dx*dx) ) * (left - 2*mid + right);
    
}


double rhs(double CFL, double left, double mid, double right)
{

    return - CFL * left + (2.0*CFL - 1) * mid - CFL * right;
}


bool solveTridiag(double* a, double* b, double* c, double* d, double* correct, int N)
{
    
    for ( int j = 2, n = N; j < n; j++ )
    {
        b[j] = b[j] - c[j-1] * a[j]/b[j-1];
        d[j] = d[j] - d[j-1] * a[j]/b[j-1];
    }
    d[N-1] = d[N-1]/b[N-1];
    
    for (int j = N-2; j > 0 ; j--)
    {
        d[j] = ( d[j] - c[j] * d[j+1] )/b[j];
    }
    
    for (int j = 1, n = N; j < n; j++)
    {
        correct[j] = d[j];
    }
    

    return true;
}