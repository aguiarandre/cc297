/**
 * helpers.c
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Implementa a execução de cada um dos 5 solvers:
 * Jacobi, GS, SOR, LGS, SLOR.
 * 
 */

 #include <stdio.h>
 #include <stdbool.h>
 #include <math.h>      /** fabs, log10 */
 #include <sys/resource.h>
 #include <sys/time.h>
 
 #include "helpers.h"
 #include "definitions.h"
 #include "solvers.h"
 
 
/**
 * Executa o ciclo de solução de Jacobi
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return Tempo de execução do ciclo de Jacobi.
 * 
 */ 
bool solveJacobi( solution* jacobi )
{
    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( jacobi ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */


    applyIC( jacobi );
    applyBC( jacobi );

    
    while ( !checkRes( jacobi, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateJacobi( jacobi );
        applyBC( jacobi );  
        
        nIterations += 1;
    }
    

    jacobi->nIterations = nIterations;
    
    return true;
    
}


/** 
 * Executa o ciclo de solução de Gauss-Seidel
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return 'true' se tudo ok, 'false' se não.
 * 
 */ 
bool solveGS( solution* gs )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( gs ) )
    {
        printf("Erro ao se inicializar a solução de Gauss Seidel.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */


    applyIC( gs );
    applyBC( gs );

    
    while ( !checkRes( gs, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateGS( gs );
        applyBC( gs );  
        
        nIterations += 1;
    }
    
    gs->nIterations = nIterations;
    
    return true;
    
}


/** TODO
 * Executa o ciclo de solução de SOR
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return 'true' se tudo ok, 'false' se não.
 * 
 */ 
bool solveSOR( solution* sor )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( sor ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    applyIC( sor );
    applyBC( sor );

    
    while ( !checkRes( sor, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateSOR( sor );
        applyBC( sor );  
        
        nIterations += 1;
    }
    
    sor->nIterations = nIterations;
    
    return true;
    
}




/**
 * Calcula o resíuo e Verifica se a solução convergiu.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis de interesse para calcular o resíduo.
 * 
 * @return 'true' se convergiu, 'false' se não.
 * 
 */ 
bool checkRes( solution * sol , int nIterations)
{
    /**
     *
     * Obter máximo valor de resíduo
     * Isto é, obter a norma infinito.
     * 
     */
    if (!sol )
    {
        printf("Erro ao enviar o ponteiro para checkRes() \n");
        return false;
    }
    
    
    /**
     * Cálculo do resíduo.
     */ 
    
    for (int i = 1, imax = IMAX-1 ; i < imax; i++ )
    {
        for (int j = 1, jmax = JMAX-1 ; j< jmax ; j++)
        {
            
            sol->res[i][j] = ( 2.0 / (sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j]) ) *
                                 ( (sol->phi[i+1][j] - sol->phi[i][j]) / ((sol->mesh->x[i+1][j] - sol->mesh->x[i][j])) -
                                   (sol->phi[i][j] - sol->phi[i-1][j]) / ((sol->mesh->x[i][j] - sol->mesh->x[i-1][j]))  ) +
                             ( 2.0 / (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1]) ) *
                                 ( (sol->phi[i][j+1] - sol->phi[i][j]) / ((sol->mesh->y[i][j+1] - sol->mesh->y[i][j] )) -
                                   (sol->phi[i][j] - sol->phi[i][j-1]) / ((sol->mesh->y[i][j] - sol->mesh->y[i][j-1] )) );
                                   
        }
    }
    
    double resMax = 0;
    
    for(int i = 1, imax = IMAX-1; i < imax ; i++)
    {
        for(int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {
            if ( fabs( sol->res[i][j] ) > resMax)
            {
                resMax = fabs( sol->res[i][j] );
            }
        }
    }
    
    /**
     * 
     * Obter a norma L2 
     * (TO DO) 
     * 
     */ 
    
    if ((nIterations%PLOT_RES == 0 ) || resMax < eps )
    {
        sol->resMax = resMax;
        printf("log10(Residuo) = %f\n", log10(resMax) );
        
    }
    
    if (resMax < eps)
    {
        printf("log10(Residuo) = %f\n", log10(resMax) );
        printf("\nConvergiu em %d iterações.\n", nIterations);
        sol->convergiu = true;
        sol->resMax = resMax;
        return true;
    }
    
    if ( !writeRes(resMax, nIterations, sol->name) )
    {
        printf("Erro ao se escrever o resíduo no arquivo.\n");
        return false;
    }
    
    return false;
}


/**
 * Aplica 1 iteração do metodo de Jacobi.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateJacobi ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Jacobi\n");
        return false;
    }
    
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
       {
        
           sol->correction[i][j] = - sol->res[i][j] / sol->mesh->N[i][j];
           
       }
   }
   
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
       {   
           sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
       }
   }
   
    return true;
}

/** TODO
 * Aplica 1 iteração do metodo de Gauss Seidel.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateGS ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    

   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       sol->correction[i][1] = ( - sol->res[i][1] 
                                 - sol->correction[i-1][1]/(sol->mesh->dxdx[i][1])  ) / ( 1 / sol->mesh->dydy[i][1] + sol->mesh->N[i][1] ) ;
       
       for (int j = 2, jmax = JMAX-1; j < jmax ; j++)
       {
          
           sol->correction[i][j] = ( - sol->res[i][j] 
                                     - sol->correction[i-1][j]/(sol->mesh->dxdx[i][j]) 
                                     - sol->correction[i][j-1]/(sol->mesh->dydy[i][j] ) ) / sol->mesh->N[i][j] ;
           
       }
   }
   
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
    

    return true;
}

bool iterateSOR ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    
    
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       sol->correction[i][1] = ( - sol->res[i][1] 
                                 - sol->correction[i-1][1]/(sol->mesh->dxdx[i][1])  ) / ( 1 / sol->mesh->dydy[i][1] + sol->mesh->N[i][1] ) ;
                                 
       for (int j = 2, jmax = JMAX-1; j < jmax ; j++)
       {
           
           sol->correction[i][j] = ( - sol->res[i][j] 
                                     - sol->correction[i-1][j]/(sol->mesh->dxdx[i][j]) 
                                     - sol->correction[i][j-1]/(sol->mesh->dydy[i][j] ) ) / sol->mesh->N[i][j] ;
           
       }
   }
   
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
    

    return true;
}

double dfdx(double x)
{
    return 2*th - 4*th*x;
}
 
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
        sol->vy[i][0] = 2*uInf*(2*th - 4*th*sol->mesh->x[i][0]) - sol->vy[i][1] ;
    }

   
    /** OVER the airfoil (mid-line) */
    double xRight,xLeft,phiTop,phiBottom, phiLeft, phiRight;

    for (int i = ILE, n = ITE+1; i < n; i++ )
    {
       
        xRight = (sol->mesh->x[i+1][0] + sol->mesh->x[i][0]) /2;
        xLeft =  (sol->mesh->x[i][0] + sol->mesh->x[i-1][0]) /2;
        
                
        phiTop =    (sol->phi[i][1] + sol->phi[i-1][1] )/2;
        phiBottom = (sol->phi[i][0] + sol->phi[i-1][0] )/2;

        phiLeft =   (phiTop + phiBottom) / 2;

        phiTop =    (sol->phi[i+1][1] + sol->phi[i][1] )/2;
        phiBottom = (sol->phi[i+1][0] + sol->phi[i][0] )/2;
        
        phiRight =  (phiTop + phiBottom) / 2;
        

        sol->velocity[0][i] = ( phiRight - phiLeft ) / ( xRight - xLeft  );
        sol->velocity[1][i] = ( sol->phi[i][1] - sol->phi[i][0] ) / (sol->mesh->y[i][1] - sol->mesh->y[i][0]);

    }
    
    /** Cp all over the place! */
    for (int i = 0; i < IMAX; i++)
    {
        for (int j = 0 ; j < JMAX; j++)
        {
            sol->cp[i][j] = 1 - ( sol->vx[i][j]*sol->vx[i][j] + sol->vy[i][j]*sol->vy[i][j] )/( uInf*uInf );    
        }
    }



} 
