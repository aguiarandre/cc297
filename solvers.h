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
bool iterateAF1( solution * ); 
bool iterateAF2( solution * ); 

bool checkRes( solution *, int );   /** Calcula o residuo e verifica se convergiu */

bool solveJacobi( solution* );    /** Executa o ciclo de solução de Jacobi */
bool solveGS( solution* );        /** Executa o ciclo de solução de Gauss Seidel */
bool solveSOR( solution* );       /** Executa o ciclo de solução de SOR */
bool solveLGS( solution* );
bool solveSLOR( solution* );
bool solveAF1( solution* );
bool solveAF2( solution* );

void calcVelocity( solution * );

#endif // SOLVERS_H 