#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define c 1.0
#define L 1.0
#define JMAX 10
#define NMAX 50

#define NBS1 1
#define NBS2 2
#define NBS3 3
#define ALPHA 0.5
#define eps 1e-8

void leapfrog ( double *  ,double *, double );
void laxWendroff ( double * u , double * uOld, double cfl);
void applyBC( int tag, double * u ,double);
bool writeSol (char * name, double * u, int n, double cfl);
bool writeNorm (char* name, int n, double G, double norm);
double amplificationFactor( double norm, double normOld);

double calculateNorm (double * u , double dx);



int main( int argc, char* argv[])
{
    system("rm Caso*");
    
    double dx = L / JMAX;
    double dt, cfl, norm, normOld, G;
    char * name;
    
    /** Mesh generation */
    
    double * x = calloc( JMAX+1 , sizeof(double));
    for ( int j = 0; j < JMAX+1; j++)
    {
        x[j] = x[0] + (double) j * dx;
    }
    
    /** Initial Perturbation */
    
    double * u = calloc( (JMAX+1) , sizeof(double) );
    double * uOld = calloc( (JMAX+1) , sizeof(double) );

    double * uInit = calloc( (JMAX+1) , sizeof(double) );
    
    uInit[JMAX] = 1.0;      
    
    
    /**
     * Caso 1.
     * Leapfrog com NBS (a)
     * cfl = 1
     */
     
    cfl = 1.0;
    /** Calculate dt */
    dt = cfl * dx / c;
    
    /** Apply IC */
    u[JMAX] = uInit[JMAX];
    uOld[JMAX] = uInit[JMAX];
    
    /** output */
    name = "Caso1";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)
    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        leapfrog( u , uOld, cfl );

        
        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

        
        
    }
   
   
   
    /**
     * Caso 2.
     * Leapfrog com NBS (a)
     * cfl = 0.5
     */
     
    cfl = 0.5;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
    for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    
    /** output */
    name = "Caso2";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)
    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        leapfrog( u , uOld, cfl );

        
        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }
        
    }   
    

    /**
     * Caso 3.
     * Leapfrog com NBS (b)
     * cfl = 0.5
     */
     
    cfl = 0.5;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso3";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        leapfrog( u , uOld, cfl );
  
      if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }



    }       
   
    /**
     * Caso 4.
     * Leapfrog com NBS (b)
     * cfl = 1.0
     */
     
    cfl = 1.0;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso4";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        leapfrog( u , uOld, cfl );
 
        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }
      
    }      
        
    /**
     * Caso 5.
     * Leapfrog com NBS (b)
     * cfl = 1.1
     */
     
    cfl = 1.1;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso5";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        leapfrog( u , uOld, cfl );

        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    /**
     * Caso 6.
     * Lax-Wendroff com NBS (c), alpha = 0 -> portanto, NBS (b) 
     * cfl = 0.5
     */
     
    cfl = 0.5;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso6";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        laxWendroff( u , uOld, cfl );

        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    /**
     * Caso 7.
     * Lax-Wendroff com NBS (c), alpha = 0 -> portanto, NBS (b) 
     * cfl = 1.0
     */
     
    cfl = 1.0;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso7";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)
    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        laxWendroff( u , uOld, cfl );

        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    
    /**
     * Caso 8.
     * Lax-Wendroff com NBS (c), alpha = 0 -> portanto, NBS (b) 
     * cfl = 1.1
     */
     
    cfl = 1.1;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso8";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        laxWendroff( u , uOld, cfl );

        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    
    
    /**
     * Caso 9.
     * Lax-Wendroff com NBS (c), alpha = 0 -> portanto, NBS (b) 
     * cfl = 1.2
     */
     
    cfl = 1.2;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso9";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        laxWendroff( u , uOld, cfl );

        if ( n == 10 || n == 30)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    
    /**
     * Caso 10.
     * Lax-Wendroff com NBS (c), alpha = 0 -> portanto, NBS (b) 
     * cfl = 0.5
     */
     
    cfl = 0.5;
    /** Calculate dt */
    dt = cfl * dx / c;
    /** Apply IC */
   for (int j = 0; j < JMAX+1; j++)
    {
        u[j] = uInit[j];
        uOld[j] = uInit[j];
    }
    /** output */
    name = "Caso10";
    
    
    for ( int n = 1, nmax = NMAX+1; n < nmax ; n++)

    {  
        if (!writeSol( name, u , n, cfl) )
        {
            printf("Erro I/O\n");
            return -1;
        }
        
        
        laxWendroff( u , uOld, cfl );

        if ( n == 10 || n == 30 || n == 50)
        {
            norm = calculateNorm(u, dx);
            normOld = calculateNorm(uOld, dx);
        
        
            if ( norm > eps && normOld > eps)
            {
                G = amplificationFactor(norm, normOld);
            }
            else
            {
                G = 0.0;
            }
            
            if (!writeNorm( name, n, G , norm) )
            {
                printf("Erro I/O\n");
                return -1;
            }
        }

    } 
    
    
     
   
    
    return 0;
}

void applyBC( int tag, double * u , double cfl)
{
    switch(tag)
    {
        case 1 :
        
        u[JMAX] = u[JMAX-1];
        
        break;
        
        case 2 :
        u[JMAX] = u[JMAX-1] - cfl * (u[JMAX] - u[JMAX-1]);
        break;
        
        case 3 :
        u[JMAX] = u[JMAX-1] - cfl * (u[JMAX] - u[JMAX-1]) - cfl*ALPHA * (u[JMAX] - 2.0*u[JMAX-1] + u[JMAX-2]) ;

        break;
    
    }
    
}

void leapfrog ( double * u , double * uOld, double cfl)
{
    double * uNew = malloc( JMAX * sizeof(double));

    for(int i = 1; i < JMAX; i++)
    {
        uNew[i] = uOld[i] - cfl * (u[i+1] - u[i-1]);
    }
    
    applyBC( NBS1 , uNew, cfl );
    
    for(int j = 1; j < JMAX; j++)
    {
        uOld[j] = u[j];
        u[j] = uNew[j];
    }
    
    free(uNew);
}

void laxWendroff ( double * u , double * uOld, double cfl)
{
    double * uNew = malloc( JMAX * sizeof(double));
    for(int i = 1; i < JMAX; i++)
    {
        uNew[i] = u[i] - (cfl/2.0) * (u[i+1] - u[i-1]) + ( ( cfl * cfl)/2.0 ) * (u[i+1] - 2.0*u[i] + u[i-1]) ;
        
    }
    
    for(int j = 1; j < JMAX; j++)
    {
        uOld[j] = u[j];
        u[j] = uNew[j];
    }
    
    free(uNew);
    
}



bool writeSol (char * name, double * u, int n, double cfl)
{
    
    FILE * fw;
    char * file = malloc ( (strlen(name) + 4)*sizeof(char*));
    strcpy(file, name);
    strcat(file , ".dat");
    
    if ( !(fw = fopen(file,"a") ) )
    {
        return false;
    }

    fprintf(fw, "%s, CFL = %.2f, n = %d\n", name, cfl, n);
    fprintf(fw, "\n");
    for (int i = 0; i<JMAX+1; i++)
        fprintf(fw, "%.10f\n", u[i]);
    fprintf(fw, "\n");    
    
    fclose(fw);
    return true;
}

double calculateNorm (double * u , double dx)
{
    double sum = 0;
    
    for (int j = 0; j < JMAX+1; j++)
    {
        sum = sum + u[j]*u[j];
    }

return sqrt (sum * dx) ;
    
}

double amplificationFactor( double norm, double normOld)
{
    return norm/normOld;
}

bool writeNorm (char* name, int n, double G, double norm)
{
   
    FILE * fw;
    char * file = malloc ( (strlen(name) + 4)*sizeof(char*));
    strcpy(file, name);
    strcat(file , ".dat");
    
    if ( !(fw = fopen(file,"a") ) )
    {
        return false;
    }

    fprintf(fw,"\nn = %d, G = %f, L2 Norm = %f\n", n, G, norm);
    fprintf(fw, "\n");  
    
    fclose(fw);
    return true;
}