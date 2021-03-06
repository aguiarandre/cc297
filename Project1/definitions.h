/**
 * helpers.h
 *
 * CC 297 - CFD
 * Projeto 1
 *
 * Decalara 'defines' comuns ao código.
 */
 
 
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define ILE 10          // Índice do 'Leading Edge'     (11o ponto em C)
#define ITE 30          // Índice do 'Trailing Edge'    (31o ponto em C)
#define IMAX 41         // Índice máximo em X           (0 até 40  em C)
#define JMAX 12         // Índice máximo em Y           (0 até 11  em C)
#define XSF 1.25        // Fator de estiramento em X    
#define YSF 1.25        // Fator de estiramento em Y


#define PI 3.14159265358979323846
                                
#define R_SOR 1.812             // Fator de sobre-relaxação de SOR
#define R_SLOR 1.88             // Fator de sobre-relaxação de SLOR

#define ALPHA 35                // Fator de aceleração de convergência (ADI)
#define OMEGA 1.95              // R_AF1    

#define ALPHA_2 3.11            // Fator de aceleração de convergência (AF-2)   
#define OMEGA_2 1.65            // R_AF2


#define PLOT_RES 5000           // De quanto em quanto tempo, plota resíduo na tela.
#define uInf 1.0                // Velocidade - Escoamento não perturbado
#define th 0.05                 // Airfoil Thickness
#define eps 1e-12                // Convergence Criterion.
#define BLOW 1e2                // Divergence Criterion.
#define MAX_ITERATIONS 20000    // Número máximo de iterações.
#define MESH_NAME "malha.dat"   // Nome do resultado a ser lido pelo Tecplot/Paraview

#endif // DEFINITIONS_H
