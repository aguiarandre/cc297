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
double solveJacobi( solution* jacobi )
{
    /** structs que guardam informações sobre tempo, do tipo 'rusage' */
    struct rusage before, after;
    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( jacobi ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return 3;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */

    getrusage(RUSAGE_SELF, &before);
    
    applyIC( jacobi );
    applyBC( jacobi );

    
    while ( !checkRes( jacobi, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateJacobi( jacobi );
        applyBC( jacobi );  
        
        nIterations += 1;
    }
    
    getrusage(RUSAGE_SELF, &after);
    
    /** Calcula o tempo de execução de Jacobi */
    return calculate(&before, &after);
    
}


/** TODO
 * Executa o ciclo de solução de Gauss-Seidel
 * 
 * @param &solution   Endereço que aponta para uma 'struct solution', que contém 
 *                    as variáveis de interesse.
 * 
 * @return Tempo de execução do ciclo de Jacobi.
 * 
 */ 
double solveGS( solution* gs )
{
    /** structs que guardam informações sobre tempo, do tipo 'rusage' */
    struct rusage before, after;
    
     /** 
     * Alocar memória às variáveis da solução.
     */
    if ( !solutionInit( gs ) )
    {
        printf("Erro ao se inicializar a solução de Jacobi.\n");
        return 3;
    }  
    
    int nIterations = 0;
    
    /** Condição inicial */
    
    applyIC( gs );
    applyBC( gs );

    getrusage(RUSAGE_SELF, &before);
    
    while ( !checkRes( gs, nIterations ) && nIterations < MAX_ITERATIONS )
    {
        iterateGS( gs );
        applyBC( gs );  
        
        nIterations += 1;
    }
    
    getrusage(RUSAGE_SELF, &after);
    
    /** Calcula o tempo de execução de GS */
    return calculate(&before, &after);
    
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
            
            sol->res[i][j] = ( 2 / (sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j]) ) *
                                 ( (sol->phi[i+1][j] - sol->phi[i][j]) / ((sol->mesh->x[i+1][j] - sol->mesh->x[i][j])) -
                                   (sol->phi[i][j] - sol->phi[i-1][j]) / ((sol->mesh->x[i][j] - sol->mesh->x[i-1][j]))  ) +
                                 ( 2/ (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1]) ) *
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
    
    if ((nIterations%100 == 0 ))
    {
        printf("log10(Residuo) = %f\n", log10(resMax) );
    }
    
    if (resMax < eps)
    {
        printf("\nConvergiu em %d iterações.\n", nIterations);
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
    
   double C, N, dx, dy;
   
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
       {
           dx = (sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j] )/2;
           dy = (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1] )/2;
           
           N = -2/(dx*dx) -2/(dy*dy);
           
           C = - sol->res[i][j] / N;
           
           sol->phi[i][j] = sol->phi[i][j] + C; 
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
    
   double C, N, dx, dy;
   
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
       {
           dx = (sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j] )/2;
           dy = (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1] )/2;
           
           N = -2/(dx*dx) -2/(dy*dy);
                       
            sol->res[i][j] = ( 2 / (sol->mesh->x[i+1][j] - sol->mesh->x[i-1][j]) ) *
                                 ( (sol->phi[i+1][j] - sol->phi[i][j]) / ((sol->mesh->x[i+1][j] - sol->mesh->x[i][j])) -
                                   (sol->phi[i][j] - sol->phi[i-1][j]) / ((sol->mesh->x[i][j] - sol->mesh->x[i-1][j]))  ) +
                                 ( 2/ (sol->mesh->y[i][j+1] - sol->mesh->y[i][j-1]) ) *
                                 ( (sol->phi[i][j+1] - sol->phi[i][j]) / ((sol->mesh->y[i][j+1] - sol->mesh->y[i][j] )) -
                                   (sol->phi[i][j] - sol->phi[i][j-1]) / ((sol->mesh->y[i][j] - sol->mesh->y[i][j-1] )) );
           
           C = - sol->res[i][j] / N;
           
           sol->phi[i][j] = sol->phi[i][j] + C; 
           
       }
   }

    return true;
}