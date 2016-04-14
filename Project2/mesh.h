/**
 * mesh.h
 *
 * CC 297 - CFD
 * Projeto 2
 *
 * Decalara a funcionalidade da malha.
 */

#ifndef MESH_H
#define MESH_H

/**
 * Definições para as variáveis de malha.
 */ 

/** Struct que armazena os nós da malha */
 
typedef struct domain
{
    double** x;  /** Real grid */  
    double** y;
    
    double ** xRef; /** Local reference grid */
    double ** yRef;
    
    double A;
    double B;
    double C;
    
    double psi;
    double phi;
    
} domain;


/**
 * Funções que auxiliam na geração da malha.
 */
 
bool parabolicMesh(domain *); 
bool innerBC(domain *);
bool outerBC(domain *);
bool referenceGrid(domain *, int ); 

bool meshCreate( domain * );
bool stepMesh( domain *, int);
double airfoilConfig (double x);
 
 
bool meshInit(domain *);        /** Aloca memória à malha (à 'struct domain') */
bool meshCreate(domain *);      /** Popula a malha com os valores definidos em "definitions.h" */
bool meshDestroy(domain *);     /** Libera o espaço de memória da 'struct domain'*/

#endif // MESH_H