/**
 * projeto1.c
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Implementa Solvers para a equação de Laplace. Resolve o problema de um 
 * aerofólio biconvexo através de uma modelagem de pequenas perturbações.
 * 
 * Mais especificações em :
 * http://www.comp.ita.br/~azevedo/material_distribuido/hand04.pdf
 * 
 */
/** Count Lines of Code find . \( -iname \*.c -o -iname \*.h \) -exec wc -l '{}' \+ */

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
           timeSOR = 0.0, timeMesh = 0.0,
           timeLGS = 0.0, timeSLOR = 0.0,
           timeAF1 = 0.0, timeAF2 = 0.0;
    
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
       
    solution jacobi = { &mesh,      /* aponta qual malha ele deve usar. */
                        NULL,       /* inicializando phi = NULL */
                        NULL,       /* inicializando res = NULL */
                        NULL,       /* inicializando correction = NULL */
                        NULL,       /* inicializando velocity = NULL */
                        NULL,       /* inicializando vx = NULL */
                        NULL,       /* inicializando vy = NULL */
                        NULL,       /* inicializando cp = NULL */
                        0,          /* inicializando nIterations = 0 */
                        0,          /* inicializando resMax = 0 */
                        "jacobi",   /* nome */ 
                        false };    /* convergiu? */
    
    if ( !solveJacobi( &jacobi ) )
    {
        printf("Erro em Jacobi.\n");
        return 3;
    }    

    calcVelocity( &jacobi );
    calcCP( &jacobi );

    getrusage(RUSAGE_SELF, &after);
    timeJacobi = calculate(&before, &after);


    /************ GAUSS SEIDEL ***************/
    getrusage(RUSAGE_SELF, &before);

     /** Inicializando solução - GAUSS SEIDEL */ 
    
    printf("\nIniciando Solução de Gauss Seidel para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    

    solution gaussSeidel = { &mesh, NULL, NULL,NULL,NULL, NULL,NULL,NULL, 0, 0, "gaussSeidel", false };
                    
    if ( !solveGS( &gaussSeidel ) )
    {
        printf("Erro em Gauss Seidel\n");
        return 3;
    }

    calcVelocity( &gaussSeidel );
    calcCP( &gaussSeidel );

    getrusage(RUSAGE_SELF, &after); 
    timeGS = calculate(&before, &after);
    
    /********** LINE GAUSS SEIDEL *************/

     getrusage(RUSAGE_SELF, &before);

     /** Inicializando solução - LINE GAUSS SEIDEL */ 
    
    printf("\nIniciando Solução de Line Gauss Seidel para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);

    solution lgs = { &mesh, NULL, NULL,NULL,NULL, NULL,NULL,NULL, 0, 0, "lgs", false };
                    
    if ( !solveLGS( &lgs ) )
    {
        printf("Erro em Line Gauss Seidel\n");
        return 3;
    }

    calcVelocity( &lgs );
    calcCP( &lgs );
    
    getrusage(RUSAGE_SELF, &after); 
    timeLGS = calculate(&before, &after);
   
    /***************** SOR *******************/
    
    getrusage(RUSAGE_SELF, &before);

    printf("\nIniciando Solução de SOR para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    solution SOR = { &mesh, NULL, NULL, NULL,NULL,NULL,NULL, NULL, 0, 0, "SOR", false };

    if( !solveSOR( &SOR ) ) 
    {
        printf("Erro em SOR\n");
        return 3;
    }

    calcVelocity( &SOR );
    calcCP( &SOR );
    
    getrusage(RUSAGE_SELF, &after);
    timeSOR = calculate(&before, &after);
    
    /************* LINE SOR ******************/
    
     getrusage(RUSAGE_SELF, &before);

    /** Inicializando solução - LINE GAUSS SEIDEL */ 
    printf("\nIniciando Solução de Line-SOR para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    

    solution SLOR = { &mesh, NULL, NULL,NULL,NULL, NULL,NULL,NULL, 0, 0, "SLOR", false };
                    
    if ( !solveSLOR( &SLOR ) )
    {
        printf("Erro em SLOR \n");
        return 3;
    }

    calcVelocity( &SLOR );
    calcCP( &SLOR );
    
    getrusage(RUSAGE_SELF, &after); 
    timeSLOR = calculate(&before, &after);
  
     /************* AF1 SCHEME ******************/
    
     getrusage(RUSAGE_SELF, &before);

    /** Inicializando solução - LINE GAUSS SEIDEL */ 
    printf("\nIniciando Solução de AF1 para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    
    solution AF1 = { &mesh, NULL, NULL,NULL,NULL, NULL,NULL,NULL, 0, 0, "AF1", false };
                    
    if ( !solveAF1( &AF1 ) )
    {
        printf("Erro em AF1 \n");
        return 3;
    }

    calcVelocity( &AF1 );
    calcCP( &AF1 );
    
    getrusage(RUSAGE_SELF, &after); 
    timeAF1 = calculate(&before, &after);
    
    
    /************* AF2 SCHEME ******************/
    
     getrusage(RUSAGE_SELF, &before);

     /** Inicializando solução - LINE GAUSS SEIDEL */ 
    
    printf("\nIniciando Solução de AF2 para th = %.0f\%% e Uinf = %.1f m/s.\n", th*100, uInf);
    

    solution AF2 = { &mesh, NULL, NULL,NULL,NULL,NULL,NULL,NULL, 0, 0, "AF2",false };
                    
    if ( !solveAF2( &AF2 ) )
    {
        printf("Erro em AF2 \n");
        return 3;
    }

    calcVelocity( &AF2 );
    calcCP( &AF2 );

    getrusage(RUSAGE_SELF, &after); 
    timeAF2 = calculate(&before, &after);
  

  
  
    /************ END PROGRAM ****************/
        
        
    printf("\n");
    printf("Mesh Gen:   Time, %.3f s\n", timeMesh);
    printf("Jacobi:     Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeJacobi, jacobi.convergiu?"sim":"não", log10(jacobi.resMax), jacobi.nIterations);
    printf("GS:         Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeGS,     gaussSeidel.convergiu?"sim":"não", log10(gaussSeidel.resMax), gaussSeidel.nIterations );
    printf("SOR:        Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeSOR,    SOR.convergiu?"sim":"não" , log10(SOR.resMax), SOR.nIterations);
    printf("LGS:        Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeLGS,    lgs.convergiu?"sim":"não",  log10(lgs.resMax), lgs.nIterations );
    printf("SLOR:       Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeSLOR,   SLOR.convergiu?"sim":"não", log10(SLOR.resMax), SLOR.nIterations );
    printf("AF1:        Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeAF1,    AF1.convergiu?"sim":"não", log10(AF1.resMax), AF1.nIterations );
    printf("AF2:        Time, %.3f s, Convergiu:  %s, resMax = %06.3f - %d Iterations\n", timeAF2 ,   AF2.convergiu?"sim":"não", log10(AF2.resMax), AF2.nIterations);
    printf("Total time: %f", timeJacobi+ timeGS + timeSOR + timeLGS+ timeSLOR+ timeAF1+ timeAF2);

    printf("\n");

    if ( !writeSolution( "cp.dat", &jacobi, &gaussSeidel, &SOR, &lgs , &SLOR, &AF1, &AF2) )
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

    if ( !solutionDestroy( &lgs ) )
    {
        printf("Erro ao se liberar a memória da solução - LGS");
        return 3;
    }
    if ( !solutionDestroy( &SLOR ) )
    {
        printf("Erro ao se liberar a memória da solução - SLOR");
        return 3;
    }    
    
    if ( !solutionDestroy( &AF1 ) )
    {
        printf("Erro ao se liberar a memória da solução - SLOR");
        return 3;
    }
    
    if ( !solutionDestroy( &AF2 ) )
    {
        printf("Erro ao se liberar a memória da solução - SLOR");
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

