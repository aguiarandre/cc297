/**
 * projeto2.c
 *
 * CC 297 - CFD
 * Projeto 2
 *
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
        printf("Uso: ./projeto2 [nomeOutput.dat] \n");
        return 1;
    }
    
    /** Nome padrão - malha.dat - ou o que o usuário escolher */
    //char* fileName = (argc == 1) ? MESH_NAME : argv[1];
 
    /** benchmarks para o cálculo do tempo*/
    double time = 0.0;
    
    domain biconvex;
    

    getrusage(RUSAGE_SELF, &before);
    
    
    meshInit( &biconvex );
    innerBC ( &biconvex );
    outerBC ( &biconvex );
    
    parabolicMesh( &biconvex );
    
    getrusage(RUSAGE_SELF, &after); 
    time = calculate(&before, &after);

    FILE * fw = fopen("geometry.dat", "w");
    
    for (int i = 0 ; i <= IMAX-1; i++)
    {
        fprintf(fw, "%f  %f,\n", biconvex.x[i][0], biconvex.y[i][0]);
    }

    fclose(fw);
    /*
    for(int i = 0; i < JMAX; i++)
  printf("x[0][%d] = %f, xRef[0][%d] = %f\n", i, biconvex.x[0][i],i, biconvex.xRef[0][i]);
  
  printf("\n");
    for(int i = 0; i < JMAX; i++)
  printf("y[23][%d] = %f, yRef[23][%d] = %f\n", i, biconvex.y[0][i],i, biconvex.yRef[0][i]);
  */
  
  writeTecplot("malha.dat", &biconvex);
  
  meshDestroy( &biconvex );
  
  
    /************ END PROGRAM ****************/

    
    /** Roda script para plotar residuos no gnuplot */    
    system("gnuplot gnuscript");
    /** That's all folks! */
    
    return 0;
}

