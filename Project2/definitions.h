/**
 * helpers.h
 *
 * CC 297 - CFD
 * Projeto 2
 *
 * Decalara 'defines' comuns ao código.
 */
 
 
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define PI 3.14159265358979323846

#define R_FARFIELD 6.5          // Raio = 6.5 cordas
#define th 0.10                 // Airfoil Thickness
#define IMAX 93                 // CSI MAX - Azimutal
#define JMAX 15                 // ETA MAX - Radial
#define CHORD 1.0               
#define MESH_STRETCHING 0.1     // 10% de estiramento sobre o aerofólio.

#define R_SOR 1.812             // Fator de sobre-relaxação de SOR
#define R_SLOR 1.88             // Fator de sobre-relaxação de SLOR
#define ALPHA 35                // Fator de aceleração de convergência (ADI)
#define OMEGA 1.95              // R_AF1    
#define ALPHA_2 3.11            // Fator de aceleração de convergência (AF-2)   
#define OMEGA_2 1.65            // R_AF2


#define PLOT_RES 5000           // De quanto em quanto tempo, plota resíduo na tela.
#define eps 1e-6                // Convergence Criterion.
#define BLOW 1e2                // Divergence Criterion.
#define MAX_ITERATIONS 20000    // Número máximo de iterações.
#define MESH_NAME "malha.dat"   // Nome do resultado a ser lido pelo Tecplot/Paraview

#define Q_ESTIRAMENTO 1.15


#define ITE 1
#define ILE 1
#define uInf 1

#endif // DEFINITIONS_H
