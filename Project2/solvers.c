/**
 * helpers.c
 *
 * CC 297 - CFD
 * Projeto 2
 *
 * 
 */
    
    
 #include <stdio.h>         /** printf */
 #include <stdlib.h>        /** malloc */
 #include <stdbool.h>       /** true e false (poderia definir TRUE = 1...) */
 #include <math.h>          /** fabs, log10 */
 #include <sys/resource.h>  /** rusage*, para tempo comp. */
 #include <sys/time.h>      /** rusage*, para tempo comp. */
 
 #include "helpers.h"
 #include "definitions.h"
 #include "solvers.h"
 

/**
 * Calcula o resíuo e Verifica se a solução convergiu.
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém
 *                    as variáveis de interesse para calcular o resíduo.
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
            
            sol->res[i][j] = 1;
            
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





