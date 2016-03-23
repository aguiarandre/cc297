/**
 * projeto1.c
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Implementa um Solver para a equação de Laplace.
 * 
 * Mais especificações em :
 * http://www.comp.ita.br/~azevedo/material_distribuido/hand04.pdf
 * 
 */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <math.h>

#undef calculate
#undef getrusage

#include "mesh.h"
#include "helpers.h"
#include "solvers.h"
#include "definitions.h"


int main(int argc, char* argv[])
{
     struct rusage before, after;

    /**
     * Alguns 'sanity checks'
     */ 
    if (argc > 2)
    {
        printf("Uso: ./projeto1 [nomeOutput.dat] \n");
        return 1;
    }
    
    /** Nome padrão - malha.dat - ou o que o usuário escolher */
    char* fileName = (argc == 1) ? MESH_NAME : argv[1];
 
    /** benchmarks para o cálculo do tempo*/
    double timeJacobi = 0.0, timeGS = 0.0,
           timeSOR = 0.0, timeMesh = 0.0;
    
    getrusage(RUSAGE_SELF, &before);

    /** Inicializando o tipo de malha! */
    domain mesh;
    
    /**
     * Alocar memória à malha.
     */ 
    if ( !meshInit( &mesh ) )
    {
        printf("Erro ao alocar memória para a malha!\n");
        return 2;
    }

    /**
     * Popular as arrays desta malha com as definições do Projeto 1.
     */    
    if ( !meshCreate( &mesh ) )
    {
        printf("Erro ao se popular a malha!\n");
        return 2;
    }

    getrusage(RUSAGE_SELF, &after);
    
    timeMesh = calculate(&before, &after);
    
    

    /*************** JACOBI ****************/
    getrusage(RUSAGE_SELF, &before);
    /** Inicializando solução - JACOBI */ 
    printf("\nIniciando Solução de Jacobi para th = %.0f%% e Uinf = %.1f m/s.\n", th*100, uInf); 
       
    solution jacobi = { &mesh, 
                        NULL, 
                        NULL, 
                        NULL,
                        NULL, 
                        NULL,
                        NULL,
                        NULL, 
                        0, 
                        0, 
                        "jacobi", 
                        false };
    //jacobi.mesh = &mesh;    /* aponta qual malha ele deve usar. */
    //jacobi.name = "jacobi";
    //jacobi.convergiu = false;
    if ( !solveJacobi( &jacobi ) )
    {
        printf("Erro em Jacobi.\n");
        return 3;
    }    

    calcVelocity( &jacobi );

    getrusage(RUSAGE_SELF, &after);
    timeJacobi = calculate(&before, &after);


    /************ GAUSS SEIDEL ***************/
    getrusage(RUSAGE_SELF, &before);

     /** Inicializando solução - GAUSS SEIDEL */ 
    
    printf("\nIniciando Solução de Gauss Seidel para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    

    solution gaussSeidel = { &mesh, 
                    NULL, 
                    NULL,
                    NULL,
                    NULL, 
                    NULL,
                    NULL,
                    NULL, 
                    0, 
                    0, 
                    "gaussSeidel", 
                    false };
                    
    if ( !solveGS( &gaussSeidel ) )
    {
        printf("Erro em Gauss Seidel\n");
        return 3;
    }

    calcVelocity( &gaussSeidel );

    getrusage(RUSAGE_SELF, &after); 
    timeGS = calculate(&before, &after);
    
    /********** LINE GAUSS SEIDEL *************/

   
   
  

   
    
    /***************** SOR *******************/
     getrusage(RUSAGE_SELF, &before);
    /** Para SOR e SLOR, o que estou chamando de N muda. 
     *  Portanto calcularei N aqui para que seja calculado 
     *  apenas uma vez.
     */
    
    for (int i = 1, imax = IMAX-1; i < imax ; i++)
    {
        for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
        {
            mesh.N[i][j] = (-2 / (R_SOR*mesh.dxdx[i][j]) - 2 / (R_SOR*mesh.dydy[i][j]));        
        }
    }

    printf("\nIniciando Solução de SOR para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
 
    solution SOR = { &mesh, 
                    NULL, 
                    NULL, 
                    NULL,
                    NULL,
                    NULL,
                    NULL, 
                    NULL, 
                    0, 
                    0, 
                    "SOR", 
                    false };

    if( !solveSOR( &SOR ) ) 
    {
        printf("Erro em SOR\n");
        return 3;
    }

    calcVelocity( &SOR );
    
    getrusage(RUSAGE_SELF, &after);
    timeSOR = calculate(&before, &after);
    /************* LINE SOR ******************/
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    /************ END PROGRAM ****************/
        
    printf("\n");
    printf("Mesh:   Time, %.3f s\n", timeMesh);
    printf("Jacobi: Time, %.3f s, %d Iterations\n", timeJacobi, jacobi.nIterations);
    printf("GS:     Time, %.3f s, %d Iterations \n", timeGS, gaussSeidel.nIterations );
    printf("SOR:    Time, %.3f s, %d Iterations \n", timeSOR, SOR.nIterations );

    printf("\n");

    if ( !writeSolution( "cp.dat", &jacobi, &gaussSeidel, &SOR) )
    {
        printf("Erro ao escrever resultados.");
        return 3;
    }    
    if ( !writeTecplot( fileName , &gaussSeidel ) )
    {
        printf("Erro ao escrever a malha no arquivo.\n");
        return 3;
    }
    
    /**
    * Liberar memória alocada para a 'struct solution' 
    * para fins de evitar memory leak 
    */    
    
    if ( !solutionDestroy( &jacobi ) )
    {
        printf("Erro ao se liberar a memória da solução - Jacobi");
        return 3;
    }

    if ( !solutionDestroy( &gaussSeidel ) )
    {
        printf("Erro ao se liberar a memória da solução - Gauss Seidel");
        return 3;
    }
    if ( !solutionDestroy( &SOR ) )
    {
        printf("Erro ao se liberar a memória da solução - SOR");
        return 3;
    }


    /**
     * Liberar espaço alocados à malha 
     * para evitar memory leak.
     */    
     
    if ( !meshDestroy( &mesh ) )
    {
        printf("Erro ao liberar a memória da malha!\n");
        return 2;
    }
    
    /** Roda script para plotar residuos no gnuplot */    
    system("gnuplot gnuscript");
    
//    system("open results.png");    
    /** That's all folks! */
    
    return 0;
}

