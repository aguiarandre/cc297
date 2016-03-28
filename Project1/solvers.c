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
 #include <stdlib.h>
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


/** 
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
    
    if (resMax > BLOW)
    {
        printf("Divergiu em %d iterações!\n", nIterations);
        sol->convergiu = false;
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
       /** for j = 1, we have a different equation 
         * because we know C(i,0) = C(i,1).
         */
       sol->correction[i][1] = ( - sol->res[i][1] 
                                 - sol->correction[i-1][1]/(sol->mesh->dxdx[i][1])  ) / ( 1.0 / sol->mesh->dydy[i][1] + sol->mesh->N[i][1] ) ;
       
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
        printf("Erro ao enviar o ponteiro solution * para a iteração de SOR \n");
        return false;
    }
    
    
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       sol->correction[i][1] = ( - sol->res[i][1] 
                                 - sol->correction[i-1][1]/(sol->mesh->dxdx[i][1])  ) / ( 1.0 / sol->mesh->dydy[i][1] + (sol->mesh->N[i][1]/R_SOR) ) ;
                                 
       for (int j = 2, jmax = JMAX-1; j < jmax ; j++)
       {
           
           sol->correction[i][j] = ( - sol->res[i][j] 
                                     - sol->correction[i-1][j]/(sol->mesh->dxdx[i][j]) 
                                     - sol->correction[i][j-1]/(sol->mesh->dydy[i][j] ) ) / (sol->mesh->N[i][j]/R_SOR) ;
           
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



/** TODO
 * Executa o ciclo de solução de SOR
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return 'true' se tudo ok, 'false' se não.
 * 
 */ 
bool solveLGS( solution* lgs )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( lgs ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    applyIC( lgs );
    applyBC( lgs );

    
    while ( !checkRes( lgs, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateLGS( lgs );
        applyBC( lgs );  
        
        nIterations += 1;
    }
    
    lgs->nIterations = nIterations;
    
    return true;
    
}

/** TODO
 * Aplica 1 iteração do metodo de Line Gauss Seidel.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateLGS ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    

    double * diagA = malloc( JMAX * sizeof(double) );
    double * diagB = malloc( JMAX * sizeof(double) );
    double * diagC = malloc( JMAX * sizeof(double) );
    double * diagD = malloc( JMAX * sizeof(double) );

    double R,B,F;
    
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
    
        for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
        {
            R = sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1];
            F = sol->mesh->y[i][j+1] - sol->mesh->y[i][j];
            B = sol->mesh->y[i][j] - sol->mesh->y[i][j-1];
            
            
            diagA[j] = 2.0/(R*B);
            diagB[j] = (-2/sol->mesh->dxdx[i][j] - 2.0/(R*F) - 2.0/(R*B));
            diagC[j] = 2.0/(R*F);
            diagD[j] = - sol->res[i][j] - (1/sol->mesh->dxdx[i][j])*sol->correction[i-1][j];
        }
        
        diagB[1] = diagB[1] + diagA[1]; /* Condição de contorno implícita! */
        
        solveTridiag(diagA, diagB, diagC, diagD, sol->correction[i], JMAX);
        
   }
   
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
    
    free(diagA);
    free(diagB);
    free(diagC);
    free(diagD);
 
    

    return true;
}


/** 
 * Executa o ciclo de solução de SOR
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return 'true' se tudo ok, 'false' se não.
 * 
 */ 
bool solveSLOR( solution* SLOR )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( SLOR ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    applyIC( SLOR );
    applyBC( SLOR );

    
    while ( !checkRes( SLOR, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateSLOR( SLOR );
        applyBC( SLOR );  
        
        nIterations += 1;
    }
    
    SLOR->nIterations = nIterations;
    
    return true;
    
}

/** TODO
 * Aplica 1 iteração do metodo de Line Gauss Seidel.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateSLOR ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    

    double * diagA = malloc( JMAX * sizeof(double) );
    double * diagB = malloc( JMAX * sizeof(double) );
    double * diagC = malloc( JMAX * sizeof(double) );
    double * diagD = malloc( JMAX * sizeof(double) );

    double R,B,F;
    
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
    
        for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
        {
            R = (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1]) * R_SLOR ;
            F = sol->mesh->y[i][j+1] - sol->mesh->y[i][j];
            B = sol->mesh->y[i][j] - sol->mesh->y[i][j-1];
            
            
            diagA[j] = 2.0/(R*B);
            diagB[j] = (-2/(R_SLOR*sol->mesh->dxdx[i][j]) - 2.0/(R*F) - 2.0/(R*B));
            diagC[j] = 2.0/(R*F);
            diagD[j] = - sol->res[i][j] - (1.0/sol->mesh->dxdx[i][j])*sol->correction[i-1][j];
        }
        
        diagB[1] = diagB[1] + diagA[1]; /* Condição de contorno implícita! */
        
        solveTridiag(diagA, diagB, diagC, diagD, sol->correction[i], JMAX);
        
   }
   
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
    
    free(diagA);
    free(diagB);
    free(diagC);
    free(diagD);
 
    

    return true;
}



/** 
 * Executa o ciclo de solução de AF1
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return 'true' se tudo ok, 'false' se não.
 * 
 */ 
bool solveAF1( solution* AF1 )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( AF1 ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    applyIC( AF1 );
    applyBC( AF1 );

    
    while ( !checkRes( AF1, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateAF1( AF1 );
        applyBC( AF1 );  
        
        nIterations += 1;
    }
    
    AF1->nIterations = nIterations;
    
    return true;
    
}

/** TODO
 * Aplica 1 iteração do metodo de Line Gauss Seidel.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateAF1 ( solution * sol )
{
    if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    

    double * diagAx = malloc( IMAX * sizeof(double) );
    double * diagBx = malloc( IMAX * sizeof(double) );
    double * diagCx = malloc( IMAX * sizeof(double) );
    double * diagDx = malloc( IMAX * sizeof(double) );
    
    double ** f = calloc( IMAX , sizeof(double*) ); 
    for ( int i = 0; i < JMAX; i++)
    {
        f[i] = calloc( IMAX, sizeof(double) );
    }
    
    double * diagAy = malloc( JMAX * sizeof(double) );
    double * diagBy = malloc( JMAX * sizeof(double) );
    double * diagCy = malloc( JMAX * sizeof(double) );
    double * diagDy = malloc( JMAX * sizeof(double) );
    
    double Rx,Bx,Fx;
    double Ry, By, Fy;
    
    
    /** 
     * 
     * PASSO 1 
     * 
     */
    

    
    for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
    {
        for (int i = 1, imax = IMAX-1; i < imax ; i++)
        {
            Rx = sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j] ;
            Fx = sol->mesh->x[i+1][j] - sol->mesh->x[i][j];
            Bx = sol->mesh->x[i][j]   - sol->mesh->x[i-1][j];
            
            
            diagAx[i] = -2.0/(Rx*Bx);
            diagBx[i] = (ALPHA + 2.0/(Rx*Fx) + 2.0/(Rx*Bx));
            diagCx[i] = -2.0/(Rx*Fx);
            diagDx[i] = ALPHA*OMEGA*sol->res[i][j] ;
        }
        
        
        solveTridiag(diagAx, diagBx, diagCx, diagDx, f[j], IMAX);
        
            
    }
  
    /**
     * 
     * PASSO 2
     * 
     */ 
  
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
    
        for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
        {
   
            Ry = sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1] ;
            Fy = sol->mesh->y[i][j+1] - sol->mesh->y[i][j];
            By = sol->mesh->y[i][j]   - sol->mesh->y[i][j-1];
            
            diagAy[j] = -2.0/(Ry*By);
            diagBy[j] = (ALPHA + 2.0/(Ry*Fy) + 2.0/(Ry*By));
            diagCy[j] = -2.0/(Ry*Fy);
            diagDy[j] = f[j][i]; 
            
        }
        
        diagBy[1] = diagBy[1] + diagAy[1]; /* Condição de contorno implícita! */

        solveTridiag(diagAy, diagBy, diagCy, diagDy, sol->correction[i], JMAX);
        
    }
    
    
    /**
     * Solution
     */ 
    
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
         
    free(diagAx);
    free(diagBx);
    free(diagCx);
    free(diagDx);
    
    free(diagAy);
    free(diagBy);
    free(diagCy);
    free(diagDy);
    
    for  (int j = 0; j < IMAX; j++)
    {
        free(f[j]);
    }
    free(f);
  
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
bool solveAF2( solution* AF2 )
{

    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( AF2 ) )
    {
        printf("Erro ao se inicializar a solução de AF2.\n");
        return false;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    applyIC( AF2 );
    applyBC( AF2 );

    
    while ( !checkRes( AF2, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateAF2( AF2 );
        applyBC( AF2 );  
        
        nIterations += 1;
    }
    
    AF2->nIterations = nIterations;
    
    return true;
    
}

/** TODO
 * Aplica 1 iteração do metodo de Line Gauss Seidel.
 * 
 * @param &solution Endereço que aponta para uma 'struct solution', que contém
 *                  as variáveis necessárias para atualização da solução.
 * 
 * @return 'true' se efetuou corretamente esta iteração, 'false' se não.
 * 
 */ 
bool iterateAF2 ( solution * sol )
{
       if (!sol)
    {
        printf("Erro ao enviar o ponteiro solution * para a iteração de Gauss Seidel\n");
        return false;
    }
    

    double ** f = calloc( IMAX , sizeof(double*) ); 
    for ( int i = 0; i < IMAX; i++)
    {
        f[i] = calloc( JMAX, sizeof(double) );
    }
    
    double * diagAy = malloc( JMAX * sizeof(double) );
    double * diagBy = malloc( JMAX * sizeof(double) );
    double * diagCy = malloc( JMAX * sizeof(double) );
    double * diagDy = malloc( JMAX * sizeof(double) );
    
    double Fx, Bx;
    double Ry, By, Fy;
    
    
    /** 
     * 
     * PASSO 1 
     * 
     */
    
    for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
    {
        for (int i = IMAX-2; i > 0; i--)
        {
            Fx = sol->mesh->x[i+1][j] - sol->mesh->x[i][j];
            
            f[i][j] = (ALPHA_2*OMEGA_2*sol->res[i][j] + (f[i+1][j])/Fx ) / (ALPHA_2 + 1/Fx);
        }

    }
    
    
    /** 
     * 
     * PASSO 2 
     * 
     */
     
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
    
        for ( int j = 1, jmax = JMAX -1; j < jmax ; j++ )
        {
   
            Bx = sol->mesh->x[i][j] - sol->mesh->x[i-1][j];
            
            Ry = sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1] ;
            Fy = sol->mesh->y[i][j+1] - sol->mesh->y[i][j];
            By = sol->mesh->y[i][j]   - sol->mesh->y[i][j-1];
            
            diagAy[j] = -2.0/(Ry*By);
            diagBy[j] = ( (ALPHA_2/Bx) + 2.0/(Ry*By) + 2.0/(Ry*Fy));
            diagCy[j] = -2.0/(Ry*Fy);
            diagDy[j] = f[i][j] + (ALPHA_2*sol->correction[i-1][j] )/Bx; 
            
        }
        
        diagBy[1] = diagBy[1] + diagAy[1]; /* Condição de contorno implícita! */

        solveTridiag(diagAy, diagBy, diagCy, diagDy, sol->correction[i], JMAX);
        
    }
    
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {   
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {       
            sol->phi[i][j] = sol->phi[i][j] + sol->correction[i][j]; 
        }
    }
         

    free(diagAy);
    free(diagBy);
    free(diagCy);
    free(diagDy);
    
    for  (int j = 0; j < IMAX; j++)
    {
        free(f[j]);
    }
    free(f);
  
    return true;
}