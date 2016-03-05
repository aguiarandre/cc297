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


#define uInf 1.0                // Velocidade - Escoamento não perturbado
#define th 0.05                 // Airfoil Thickness
#define eps 1e-14               // Convergence Criterion.
#define MAX_ITERATIONS 30000    // Número máximo de iterações.
#define MESH_NAME "malha.dat"   // Nome do resultado a ser lido pelo Tecplot/Paraview

#endif // DEFINITIONS_H