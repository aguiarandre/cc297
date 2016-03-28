/**
 * mesh.c
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Implementa a funcionalidade da malha.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
 #include <sys/resource.h>
 #include <sys/time.h>

#include "mesh.h"
#include "definitions.h"




/**
 * Aloca a memória para os termos de malha na heap.
 * 
 * @param &mesh     Endereço que aponta para uma 'struct domain', que contém 
 *                  ambas as arrays de arrays x e y.
 * 
 * @return 'true' se alocou corretamente x e y; 'false' se não.
 * 
 */ 

bool meshInit( domain* mesh )
{
    if( !mesh )
    {
        printf("Algo está errado com o ponteiro passado para meshInit() \n ");
        return false;
    }    
    // Alocar memória em x
    if (! (mesh->x = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->y = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->dxdx = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->dydy = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }    
    if (! (mesh->N = malloc( IMAX * sizeof(double*)) ) )
    {
        return false;
    }

    // Alocar memória em y.
    for (int i = 0; i < IMAX ; i++)
    {
        if (! (mesh->x[i] = malloc ( JMAX * sizeof(double))  ) ) 
        {
            return false;
        }
        if (! (mesh->y[i] = malloc ( JMAX * sizeof(double)) ) )
        { 
            return false;
        }
        if (! (mesh->dxdx[i] = malloc ( JMAX * sizeof(double)) ) )
        { 
            return false;
        }
        if (! (mesh->dydy[i] = malloc ( JMAX * sizeof(double)) ) )
        { 
            return false;
        }
        if (! (mesh->N[i] = malloc ( JMAX * sizeof(double)) ) )
        { 
            return false;
        }        
        
    }
    
    return true;
}

/**
 * Popula os elementos desta malha com valores.
 * 
 * @param &mesh     Endereço que aponta para uma 'struct domain', que contém 
 *                  ambas as arrays de arrays x e y.
 * 
 * @return 'true' se populou corretamente, 'false' se não.
 */ 

bool meshCreate( domain* mesh)
{
    if( !mesh )
    {
        printf("Algo está errado com o ponteiro passado para meshCreate() \n ");
        return false;
    }
 
    double deltaX = 1.0 / (ITE - ILE);
    
    /** Loop para Y */ 
    
    for (int i = 0; i < IMAX; i++)
    {
        mesh->y[i][0] = ( -1 ) * deltaX/2;
        mesh->y[i][1] = (  1 ) * deltaX/2;
        
        for (int j = 2; j < JMAX ; j++)    
        {
            mesh->y[i][j] = mesh->y[i][j-1] + ( mesh->y[i][j-1] - mesh->y[i][j-2] )*YSF;
        }
    }
 
    /**
     * Loops para X 
     * Nota: cuidado com numeração (começa em 0 em C, começa em 1 em Fortran)
     * O Azevedo quer que o aerofólio comece no 11o ponto e termine no 31o.
     * Por isso subtraí 1 de ITE e ILE, em 'definitions.h'.
     * 
     */
    
    /**
     * Aerofólio
     */ 
 
    for (int i = ILE , n = ITE + 1 ; i < n ; i++)
    {
        for (int j = 0; j < JMAX ; j++)
        {
            
            mesh->x[i][j] = ( i - ILE )*deltaX ;
        
        }
    }
    
     /**
     * Região Pré-aerofólio.
     */ 
 
    for (int i = ILE - 1; i > -1 ; i--)
    {
        for (int j = 0; j < JMAX ; j++)
        {
            
            mesh->x[i][j] = mesh->x[i+1][j] + ( mesh->x[i+1][j] - mesh->x[i+2][j] )*XSF;
        
        }
    }
    
     /**
     * Região Pós-aerofólio.
     */ 
 
    for (int i = ITE + 1; i < IMAX ; i++)
    {
        for (int j = 0; j < JMAX ; j++)
        {
            
            mesh->x[i][j] = mesh->x[i-1][j] + ( mesh->x[i-1][j] - mesh->x[i-2][j] )*XSF;
        
        }
    }
    
 printf("Malha criada com sucesso...\n");
 
   for (int i = 1, imax = IMAX-1; i < imax ; i++)
   {
       for (int j = 1, jmax = JMAX-1; j < jmax ; j++)
       {
            mesh->dxdx[i][j] = (mesh->x[i+1][j] - mesh->x[i-1][j] )/2 * 
                               (mesh->x[i+1][j] - mesh->x[i-1][j] )/2;
            mesh->dydy[i][j] = (mesh->y[i][j+1] - mesh->y[i][j-1] )/2 *
                               (mesh->y[i][j+1] - mesh->y[i][j-1] )/2;
            
            mesh->N[i][j] = (-2 / mesh->dxdx[i][j] - 2 / mesh->dydy[i][j]);        
       }
   }
   
   
   
 return true;   
    
}

/**
 * Libera os elementos desta malha através do comando 'free'.
 * 
 * @param &mesh     Endereço que aponta para uma 'struct domain', que contém 
 *                  ambas as arrays de arrays x e y.
 * 
 * @return 'true' se liberou corretamente, 'false' se não.
 */ 
bool meshDestroy( domain * mesh )
{
    if( !mesh )
    {
        printf("Algo está errado com o ponteiro passado para meshDestroy() \n ");
        return false;
    }
    
    for (int i = 0; i < IMAX ; i++)
    {
        free( mesh->x[i] );
        free( mesh->y[i] );
        free( mesh->dxdx[i] );
        free( mesh->dydy[i] );
        free (mesh->N[i] );
    }
    
    free(mesh->x);
    free(mesh->y);
    free(mesh->dxdx);
    free(mesh->dydy);
    free(mesh->N);
    
    return true;
}