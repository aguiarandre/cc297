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

#undef calculate
#undef getrusage

#include "mesh.h"
#include "helpers.h"
#include "solvers.h"
#include "definitions.h"


int main(int argc, char* argv[])
{

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
    double timeJacobi = 0.0, timeGS = 0.0;
    
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
    
    /*************** JACOBI ****************/
    
    /** Inicializando solução - JACOBI */ 
     
    solution jacobi;
    jacobi.mesh = &mesh;    /* aponta qual malha ele deve usar. */
    jacobi.name = "jacobi";
    timeJacobi = solveJacobi( &jacobi );
    
    /************ GAUSS SEIDEL ***************/
    
     /** Inicializando solução - GAUSS SEIDEL */ 
     
    solution gaussSeidel;
    gaussSeidel.mesh = &mesh;    /* aponta qual malha ele deve usar. */
    gaussSeidel.name = "gaussSeidel";
    timeGS = solveGS( &gaussSeidel );
    
    
    /************ END PROGRAM ****************/
    

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

    if ( !writeTecplot( fileName , &jacobi ) )
    {
        printf("Erro ao escrever a malha no arquivo.");
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
    
    /** That's all folks! */
    
    return 0;
}

