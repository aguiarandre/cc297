/**
 * mesh.h
 *
 * CC 297 - CFD
 * Projeto 1
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
    double** x;     /** Array de arrays genérica. No caso de uma 'rectilinear grid', 
                     *  poderiam ser apenas vetores */
    double** y;
    
    double** dxdx;
    double** dydy;
    
    double** N;
    
} domain;

/**
 * Funções que auxiliam na geração da malha.
 */
 
bool meshInit(domain *);        /** Aloca memória à malha (à 'struct domain') */
bool meshCreate(domain *);      /** Popula a malha com os valores definidos em "definitions.h" */
bool meshDestroy(domain *);     /** Libera o espaço de memória da 'struct domain'*/

#endif // MESH_H