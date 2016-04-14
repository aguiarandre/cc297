/**
 * helpers.c
 *
 * CC 297 - CFD
 * Projeto 2
 *
 * Implementa a funcionalidades comuns ao código.
 */
 
 
 #include <stdlib.h>
 #include <stdio.h>
 #include <stdbool.h>
 #include <math.h>
 #include <string.h>
 #include <sys/resource.h>
 #include <sys/time.h>
 
 #include "helpers.h"
 #include "definitions.h"
 
 
 
 

/**
 * Aloca memória para os elementos de solução.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis a serem (posteriormente) inicializadas.
 * 
 * @return 'true' se alocou a memória corretamente, 'false' se não.
 * 
 */ 

bool solutionInit( solution* sol )
{
    if( !sol )
    {
        printf("Algo está errado com o ponteiro passado para solutionInit() \n ");
        return false;
    }    
    
    /** Alocar memória em x para phi */
    
    if ( !(sol->phi = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }
    
    /**
     * Alocar memória em x para res - lembrar que resíduo só usa as células centrais.
     *
     */
     
    if ( !(sol->res = calloc( IMAX , sizeof(double*)) ) )
    {
        return false;
    }
    
    if ( !(sol->correction = calloc( IMAX , sizeof(double*)) ) )
    {
        return false;
    }
    
    if ( !(sol->cp = malloc ( IMAX * sizeof(double*)) ) )
    {
        return false;
    }

    if ( !(sol->velocity = malloc( 2 * sizeof(double*)) ) )
    {
        return false;
    }
    
    if ( !(sol->vx = malloc ( IMAX * sizeof(double*)) ) )
    {
        return false;
    }

    if ( !(sol->vy = malloc ( IMAX * sizeof(double*)) ) )
    {
        return false;
    }
    
    /** Alocar memória em y para phi */
    
    for (int i = 0; i < IMAX ; i++)
    {
        if ( !(sol->phi[i] = malloc ( JMAX * sizeof(double))  ) ) 
        {
            return false;
        }
        
    }

    /** Alocar memória em y para res */

    for (int i = 0; i < IMAX ; i++)
    {
        if (! (sol->res[i] = calloc ( JMAX , sizeof(double))  ) ) 
        {
            return false;
        }
        
    }

    for ( int i = 0 ; i < IMAX; i++)
    {
        if (! (sol->cp[i] = malloc( JMAX * sizeof(double)) ) )
        {
            return false;
        }
    }

    for (int i = 0; i < 2 ; i++)
    {
        if (! (sol->velocity[i] = malloc ( (IMAX) * sizeof(double))  ) ) 
        {
            return false;
        }
        
    }
    
    for (int i = 0 ; i < IMAX; i ++)
    {
        if ( !(sol->vx[i] = malloc( JMAX * sizeof(double)) ) )
        {
            return false;
        }
    }
    
    for (int i = 0 ; i < IMAX; i ++)
    {
        if ( !(sol->vy[i] = malloc( JMAX * sizeof(double)) ) )
        {
            return false;
        }
    }
    
    /** Alocar memória em y para correction */    
    /** Note que estou usando CALLOC para já impor a 'condição de contorno'
     *  na correção. Isto porque a correção em (i = 0, j) é nula pois o 
     *  resultado na condição de contorno já está correto.
     */ 
    
    for (int i = 0; i < IMAX ; i++)
    {
        if ( !(sol->correction[i] = calloc ( JMAX , sizeof(double))  ) ) 
        {
            return false;
        }
        
    }
    
    return true;
}


/**
 * Libera a memória para os elementos de solução.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis cujas memórias serão liberadas.
 * 
 * @return 'true' se liberou a memória corretamente, 'false' se não.
 * 
 */ 

bool solutionDestroy( solution* sol )
{
    if( !sol )
    {
        printf("Algo está errado com o ponteiro passado para solutionDestroy() \n ");
        return false;
    }
    
    for (int i = 0; i < IMAX ; i++)
    {
        free( sol->phi[i] );
        free( sol->res[i] );
        free( sol->correction[i] );
        free( sol->cp[i] );
        free( sol->vx[i] );
        free( sol->vy[i] );
    }
   
   free( sol->velocity[0] );
   free( sol->velocity[1] );

    
    free(sol->phi);
    free(sol->res);
    free(sol->correction);
    free(sol->cp);
    free(sol->velocity);
    free(sol->vx);
    free(sol->vy);

    return true;
}



/**
 * Aplica condições iniciais baseado nas especificações do projeto 1.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  phi - variável de interesse para se aplicar IC.
 * 
 * @return 'true' se aplicou corretamente, 'false' se não.
 * 
 */ 
bool applyIC( solution * sol )
{
    if ( !sol )
    {
        printf("Erro ao se aplicar condição inicial.\n");
        return false;
    }
    
    for (int i = 0; i < IMAX ; i++ )
    {
        for (int j = 0; j < JMAX ; j++)
        {
            sol->phi[i][j] = uInf * sol->mesh->x[i][j];
        }
    }

    return true;
}


/**
 * Aplica condições de contorno baseado nas especificações do projeto 1.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  phi - variável de interesse para se aplicar BC.
 * 
 * @return 'true' se aplicou corretamente, 'false' se não.
 * 
 */ 
bool applyBC( solution * sol )
{

    if (!sol )
    {
        printf("Erro ao se aplicar a condição de contorno!\n");
        return false;
    }

    /** BOTTOM Boundary */
    
    for (int i = 0; i < IMAX ; i++ )
    {
            sol->phi[i][0] = sol->phi[i][1];
            
    }
   
    for (int i = ILE, imax = ITE + 1; i < imax; i++)
    {
        sol->phi[i][0] = sol->phi[i][1] - ( sol->mesh->y[i][1] - sol->mesh->y[i][0] ) * uInf * ( dydx(sol->mesh->x[i][0]) );
    }

    return true;
    
        
}

/**
 * Escreve um arquivo de saída para o Tecplot (Paraview também aceita!)
 * 
 * @param fileName  'string' com o nome do arquivo de saída.
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para output da solução e malha.
 * 
 * @return 'true' se escreveu corretamente, 'false' se não.
 * 
 */ 
bool writeTecplot(char* fileName, solution * sol)
{
     FILE * fw = fopen(fileName, "w");
     if (!fw)
     {
         printf("Não foi possível abrir o arquivo.");
         return false;
     }
     
     /** TECPLOT FORMAT */
     
     fprintf(fw, "VARIABLES = \"X\", \"Y\", \"Ux\", \"Uy\", \"Cp\" \n");
     fprintf(fw, "ZONE I=%d, J=%d, F=POINT\n", JMAX, IMAX);
     
     for (int i = 0; i < IMAX ; i++)
     {
         for(int j = 0; j < JMAX ; j++)
         {
             fprintf(fw, "%f %f %f %f %f\n", sol->mesh->x[i][j], sol->mesh->y[i][j], sol->vx[i][j], sol->vy[i][j],
             sol->cp[i][j]);
         }
     }
    
     fclose(fw);
    
     FILE * fw2= fopen("velocity.dat","w");
     
     if (!fw2)
     {
         printf("Não foi possível abrir o arquivo.");
         return false;
     }
     
     for (int i = 1, imax = IMAX-1 ; i < imax ; i++) 
     {
         fprintf(fw2,"%f %f %f\n", sol->mesh->x[i][1], sol->velocity[0][i], sol->velocity[1][i]);
     }
     fclose(fw2);



     return true;
     
}


/**
 * Escreve um arquivo de saída de Resíduos
 * 
 * @param resMax        Resíduo máximo obtido em checkRes().
 * @param nIterations   Número da iteração que estou gravando.
 * @param fileName      'string' com o nome do arquivo de saída.
 * 
 * @return 'true' se escreveu corretamente, 'false' se não.
 * 
 */ 
bool writeRes( double resMax, int nIterations, char * fileName )
{
    char * file = malloc ( (strlen(fileName) + 4)*sizeof(char*));
    strcpy(file, fileName);
    strcat(file , ".dat");
    
    FILE * fw;
    if(nIterations == 0)
    {
        if (!(fw = fopen(file,"w") ) )
        {
            printf("Erro ao abrir o arquivo de Resíduo %s\n", file);
            return false;
        }
    }
    else
    {
        if (!(fw = fopen(file,"a") ) )
        {
            printf("Erro ao abrir o arquivo de Resíduo %s\n", file);
            return false;
        }
    }
    
    fprintf(fw,"%d %f\n", nIterations+1, log10(resMax));
    fclose(fw);
    
    free(file);
    
    return true;
}

/**
 * Calcula o tempo entre dois benchmarks.
 * 
 * @param rusage * b        Tempo final.
 * @param rusage * a        Tempo inicial.
 * 
 * @return Tempo total computacional.
 * 
 */ 
double calculate(const struct rusage* b, const struct rusage* a)
{
    if (b == NULL || a == NULL)
    {
        return 0.0;
    }
    else
    {
        return (    
                    (   
                        (  (a->ru_utime.tv_sec * 1000000 + a->ru_utime.tv_usec) -
                           (b->ru_utime.tv_sec * 1000000 + b->ru_utime.tv_usec)   
                        ) +
                        (   (a->ru_stime.tv_sec * 1000000 + a->ru_stime.tv_usec) -
                            (b->ru_stime.tv_sec * 1000000 + b->ru_stime.tv_usec)   
                        )   
                    )   / 1000000.0
                );
    }
}

/**
 * Escreve um arquivo de saída de Cp.
 * 
 * @param char * fileName       Nome do arquivo de saída.
 * 7 @params solution *           Solução que se deseja escrever o Cp.
 * 
 * @return 'true' se escreveu corretamente, 'false' se não.
 * 
 */ 

bool writeSolution( char* fileName, solution* jacobi, solution * gs, solution* sor, solution* lgs, solution * SLOR, solution *AF1 , solution *AF2)
{
    if (!jacobi || !gs || !sor || !lgs || !SLOR || !AF1 || !AF2)
    {
        return false;
    }
    
    FILE * fw;
    
    if (!(fw = fopen(fileName,"w") ) )
    {
        printf("Erro ao abrir o arquivo de Resultados %s\n", fileName);
        return false;
    }

    double cpJacobi, cpGS, cpSOR, cpLGS, cpSLOR, cpAF1, cpAF2, cpExact,x;
    /** Cp OVER airfoil */ 
    for (int i = ILE, n = ITE + 1; i < n ; i++)
    {
        cpJacobi = (jacobi->velocity[0][i]*jacobi->velocity[0][i] + jacobi->velocity[1][i]*jacobi->velocity[1][i])/(uInf*uInf)  -1.0;
        
        cpGS =      (gs->velocity[0][i]*gs->velocity[0][i] + gs->velocity[1][i]*gs->velocity[1][i])/(uInf*uInf)  -1.0;

        cpSOR  = (sor->velocity[0][i]*sor->velocity[0][i] + sor->velocity[1][i]*sor->velocity[1][i])/(uInf*uInf)  -1.0;
        
        cpLGS = (lgs->velocity[0][i]*lgs->velocity[0][i] + lgs->velocity[1][i]*lgs->velocity[1][i])/(uInf*uInf)  -1.0;
        
        cpSLOR = (SLOR->velocity[0][i]*SLOR->velocity[0][i] + SLOR->velocity[1][i]*SLOR->velocity[1][i])/(uInf*uInf)  -1.0;
        
        cpAF1 = (AF1->velocity[0][i]*AF1->velocity[0][i] + AF1->velocity[1][i]*AF1->velocity[1][i])/(uInf*uInf)  -1.0;
        
        cpAF2 = (AF2->velocity[0][i]*AF2->velocity[0][i] + AF2->velocity[1][i]*AF2->velocity[1][i])/(uInf*uInf)  -1.0;
        
        x = jacobi->mesh->x[i][1]*2;
        
        cpExact = - (1.0 -  (1.0 + (2.0/PI)*th*(2.0 - (x-1)*log((1.0+(x-1))/(1.0-(x-1)+0.001)) ) ) * 
        (1.0 + (2.0/PI)*th*(2.0 - (x-1)*log((1.0+(x-1))/(1.0-(x-1)+0.0001)) ) )   )  ;
        // th*th*( (3.0/(PI*PI))*(2.0 - (x+1.0)*log( (1.0+(x+1.0))/(1.0-(x+1.0)+0.001) ))*(2.0 - (x+1.0)*log( (1.0+(x+1.0))/(1.0-(x+1.0)+0.001) )) 
        //- (1.0/(PI*PI))*(  log( (1.0 +(x+1.0))/(1.0 - (x+1.0)+0.001)   ) *  log( (1.0 +(x+1.0))/(1.0 - (x+1.0)+0.001) ) - (1.0 - (x+1.0)*(x+1.0)) ) )
        
        fprintf(fw, "%f      %f       %f      %f        %f      %f      %f      %f      %f\n", 
        jacobi->mesh->x[i][0], cpJacobi, cpGS, cpSOR, cpLGS, cpSLOR, cpAF1, cpAF2, cpExact );

    }
   

    fclose(fw);
    
    return true;
}

/**
 * Calcula a derivada do perfil em questão.
 * 
 * @param x        valor da coordenada x
 * 
 * @return valor da derivada da função
 * 
 */ 
double dydx( double x)
{
    
    return 2.0*th - 4.0*th*x;
   
//return th; // rampa
// y = ( x/3 ) * ( 1 - 1/9 * x^2/d^2 + 1/54 * x^3/d^3  ) //AGARD-C   
//  return 1.0/3.0 - (1.0/(9.0*3.0))*3.0 * (x*x)/(d*d) + (1.0/(54.0*3.0))*4.0*(x*x*x)/(d*d*d); //AGARD-C
//    y= +- 0.6*[0.2969*sqrt(x) - 0.1260*x - 0.3516*x2 + 0.2843*x3 - 0.1015*x4] //naca0012
//return .6*((0.2969*.5)/sqrt(x) - 0.1260 - 2*0.3516*x + 3*0.2843*x*x - 4* 0.1015*x*x*x ); 


}

/**
 * Resolve um sistema tridiagonal
 * 
 * @param *a            Lower diagonal.
 * @param *b            Main diagonal.
 * @param *c            Upper diagonal.
 * @param *d            RHS
 * @param *correct      vector where the solution needs to be written
 * @param *N            size of the problema
 * 
 * @return 'true' se calculou corretamente, 'false' se não.
 * 
 */ 

bool solveTridiag(double* a, double* b, double* c, double* d, double* correct, const int N)
{
    
    for ( int j = 2, n = N-1; j < n; j++ )
    {
        b[j] = b[j] - c[j-1] * a[j]/b[j-1];
        d[j] = d[j] - d[j-1] * a[j]/b[j-1];
    }
    d[N-2] = d[N-2]/b[N-2];
    
    for (int j = N-3; j > 0 ; j--)
    {
        d[j] = ( d[j] - c[j] * d[j+1] )/b[j];
    }
    
    for (int j = 1, n = N-1; j < n; j++)
    {
        correct[j] = d[j];
    }
    

    return true;
}

/** 
 * Calcula a velocidade no campo inteiro e sobre o perfil.
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return NULL (Escreve na própria struct)
 * 
 */ 
 
void calcVelocity( solution * sol )
{
    /**
     * Calculate Velocity ( with BC's )
     * and then calculate Cp field.
     */

    /** Velocity Top */
    for ( int i = 0 ; i < IMAX; i ++)
    {
        sol->vx[i][JMAX-1] = uInf;
        sol->vy[i][JMAX-1] = 0.0;
    }

    for (int j = 1 ; j < JMAX ; j++ )
    {
        /** Velocity Left */
        sol->vx[0][j] = uInf;
        sol->vy[0][j] = 0.0;
        
        /** Velocity Right */
        sol->vx[IMAX-1][j] = uInf;
        sol->vy[IMAX-1][j] = 0.0;
    }

    /** Middle points */
    for (int i = 1, imax = IMAX - 1; i < imax ; i++) 
    {
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {
            sol->vx[i][j] =  ( sol->phi[i+1][j] - sol->phi[i-1][j] ) / ( ( sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j]) );
            sol->vy[i][j] =  ( sol->phi[i][j+1] - sol->phi[i][j-1] ) / ( ( sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1]) );
        }
    }

     /** Bottom */
    for (int i = 0, imax = IMAX; i < imax ; i++) 
    {
       sol->vx[i][0] = sol->vx[i][1];
       sol->vy[i][0] = - sol->vy[i][1];
    }
    
    /** When over the airfoil airfoil, the BC is ... */
    for (int i = ILE, imax = ITE+1; i< imax; i++)
    {
        sol->vy[i][0] = 2.0*uInf*(dydx( sol->mesh->x[i][0])) - sol->vy[i][1] ;
    }

   
    /** OVER the airfoil (mid-line) */
    double xRight,xLeft,phiTop,phiBottom, phiLeft, phiRight;

    for (int i = ILE, n = ITE+1; i < n; i++ )
    {
       
        /**
         * x- velocity got a little bit complicated, but ...
         * I tried to make x-average the to be the same as the 
         * y-average - that is, with one spacing only.
         */
         
        xRight = (sol->mesh->x[i+1][0] + sol->mesh->x[i][0]) /2.0;
        xLeft =  (sol->mesh->x[i][0] + sol->mesh->x[i-1][0]) /2.0;
        
                
        phiTop =    (sol->phi[i][1] + sol->phi[i-1][1] )/2.0;
        phiBottom = (sol->phi[i][0] + sol->phi[i-1][0] )/2.0;

        phiLeft =   (phiTop + phiBottom) / 2.0;

        phiTop =    (sol->phi[i+1][1] + sol->phi[i][1] )/2.0;
        phiBottom = (sol->phi[i+1][0] + sol->phi[i][0] )/2.0;
        
        phiRight =  (phiTop + phiBottom) / 2.0;
        

        sol->velocity[0][i] = ( phiRight - phiLeft ) / ( xRight - xLeft  );
        
        /** 
         * The y-velocity is as simple as the boundary condition! 
         * (because IT IS the boundary condition!)
         */ 
        sol->velocity[1][i] = ( sol->phi[i][1] - sol->phi[i][0] ) / (sol->mesh->y[i][1] - sol->mesh->y[i][0]);

    }

} 

/** 
 * Calcula o CP baseado nos valores de vx e vy de uma solução.
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return NULL (escreve na própria struct solution)
 * 
 */ 

void calcCP( solution * sol)
{
    /** Cp all over the place! */
    for (int i = 0; i < IMAX; i++)
    {
        for (int j = 0 ; j < JMAX; j++)
        {
            sol->cp[i][j] = 1 - ( sol->vx[i][j]*sol->vx[i][j] + sol->vy[i][j]*sol->vy[i][j] )/( uInf*uInf );    
        }
    }

}

bool solvePeriodicTridiag(double* a, double* b, double* c, double* d, double** r, const int N, int J)
{
    int nMax = N-1;
    
    double * n = calloc ( N , sizeof(double));
    double * m = calloc ( N , sizeof(double));

    /**
     * Solves Au = d, where A has periodic BC. 
     *  _                               _    _     _     _     _
     * |   b       c 0 0  ... 0 upperBc  |  | u(0)  |   |  d(0)  |
     * |   a       b c 0  ... 0    0     |  | u(1)  |   |  d(1)  |
     * |   0       a b c  ... 0    0     |  | u(2)  |   |  d(2)  |
     * |   0       0 a b  ... 0    0     |  |   .   |   |   .    |
     * |   .             .               |x |   .   | = |   .    |
     * |   .               .             |  |   .   |   |   .    |
     * |   .                 .           |  |       |   |        |
     * |                                 |  |       |   |        |
     * |lowerBc   0 0 0  ... 0     a     |  | u(n)  |   |  d(n)  |
     *  
     */
    
    n[0] = a[0];
    m[0] = c[nMax]; 
     
     
     
    for ( int j = 1 ; j <= nMax; j++ )
    {
        b[j] = b[j] - c[j-1] * a[j]/b[j-1];
        d[j] = d[j] - d[j-1] * a[j]/b[j-1];
        
        n[j] = 0.0  - n[j-1] * a[j]/b[j-1];
        m[j] = 0.0  - c[j-1] * m[j-1]/b[j-1];
        
        b[nMax] = b[nMax] - n[j-1] * m[j-1]/b[j-1];
        d[nMax] = d[nMax] - d[j-1] * m[j-1]/b[j-1];
        
    }
    
    d[nMax] = d[nMax]/b[nMax];
    
    for (int j = nMax-1; j >= 0 ; j--)
    {
        d[j] = ( d[j] - c[j] * d[j+1] - n[j] * d[N] )/b[j];
    }
    
    for (int j = 0; j <= nMax; j++)
    {
        r[j][J] = d[j];
    }

    free(n);
    free(m);
    
    return true;
}
