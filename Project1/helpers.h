/**
 * helpers.h
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Decalara funcionalidades comuns ao código.
 */


#ifndef HELPERS_H
#define HELPERS_H

#include "mesh.h"

/**
 * Declarações da struct que armazena a solução de cada metodo.
 */
 
typedef struct solution
{
    domain * mesh;          /** Note que este elemento apenas indica 
                             *  onde encontrar a malha a ser usada, 
                             *  não ocupando muito espaço na memória */
    double** phi;           
    double** res;
    double** correction;
    double** velocity;
    double** vx;
    double** vy;
    double** cp;

    int nIterations;
    double resMax;
    
    char* name;
    bool convergiu;
    
    
} solution;

bool solutionInit(solution* );      /** Aloca memória à solução. */
bool solutionDestroy(solution *);   /** Libera a memória da struct solution */

bool applyBC(solution *);           /** Aplica condições de contorno à solution->phi */
bool applyIC( solution *);          /** Aplica a condição inicial à solution->phi */

bool writeTecplot(char *, solution*);   /** Escreve arquivo output em formato Tecplot */
bool writeRes( double, int , char*);    /** Escreve arquivo .dat com nIterações vs Resíduo */
bool writeSolution( char*, solution*, solution*, solution*, solution* ,solution * , solution*, solution*); /** Escreve a solução em um arquivo */

double dydx(double);    /** Retorna o valor da derivada de uma função no ponto */

double calculate(const struct rusage*, const struct rusage*); /** Calcula o tempo entre dois rusage (before, after) */
bool solveTridiag(double* a, double* b, double* c, double* d, double* correct, const int);  /** Resolve um sistema tridiagonal */   

void calcVelocity( solution * );
void calcCP( solution * );

#endif // HELPERS_H 
