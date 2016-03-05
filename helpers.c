/**
 * helpers.c
 *
 * CC 297 - CFD
 * Projeto 1
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
    }
    for (int i = 0; i < IMAX ; i++)
    {
        free( sol->res[i] );
    }
    
    
    free(sol->phi);
    free(sol->res);
    
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

    double foo, dy;
    
    /** TOP Boundary */
    
    for (int i = 0, jmax = JMAX-1; i < IMAX ; i++ )
    {
        sol->phi[i][jmax] = uInf * sol->mesh->x[i][jmax];
    }
    
    /** INLET Boundary */
    
    for (int j = 0; j < JMAX ; j++)
    {
        sol->phi[0][j] = uInf * sol->mesh->x[0][j];
    }
    
    
    /** EXIT Boundary */
    

    for (int j = 0, imax = IMAX -1; j < JMAX ; j++)
    {
        sol->phi[imax][j] = uInf * sol->mesh->x[imax][j];
    }

    
    /** BOTTOM Boundary */
    
    for (int i = 0; i < IMAX ; i++ )
    {
        if (i < ILE || i > ITE)
        {
            sol->phi[i][0] = sol->phi[i][1];
        }
        else
        {
            foo = uInf * ( 2*th - 4*th*sol->mesh->x[i][0] );
            dy = sol->mesh->y[i][2] - sol->mesh->y[i][1];
            
            sol->phi[i][0] = sol->phi[i][1] - (dy) * foo;
        }
            
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
     
     fprintf(fw, "VARIABLES = \"X\", \"Y\" \n");
     fprintf(fw, "ZONE I=%d, J=%d, F=POINT\n", JMAX, IMAX);
     
     for (int i = 0; i < IMAX ; i++)
     {
         for(int j = 0; j < JMAX ; j++)
         {
             fprintf(fw, "%f %f\n", sol->mesh->x[i][j], sol->mesh->y[i][j]);
         }
     }
    
     fclose(fw);
     
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
 * Retorna tempo entre b and a.
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