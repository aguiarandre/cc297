/**
 * helpers.h
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Decalara a funcionalidade dos solvers do código.
 */
 
 
#ifndef SOLVERS_H
#define SOLVERS_H

#include "mesh.h"
#include "helpers.h"


bool iterateJacobi( solution * );   /** Executa 1 iteração do metodo de Jacobi */
bool iterateGS( solution * );       /** Executa 1 iteração do metodo de Gauss Seidel */
bool iterateSOR( solution * );      /** Executa 1 iteração do metodo de SOR  */
bool iterateLGS( solution * );      /** Executa 1 iteração do metodo de Line Gauss Seidel */
bool iterateSLOR( solution * );     /** Executa 1 iteração do metodo de SLOR  */

bool checkRes( solution *, int );   /** Calcula o residuo e verifica se convergiu */

double solveJacobi( solution* );    /** Executa o ciclo de solução de Jacobi */
double solveGS( solution* );        /** Executa o ciclo de solução de Gauss Seidel */


#endif // SOLVERS_H 