/**
 * mesh.c
 *
 * CC 297 - CFD
 * Projeto 2
 *
 * Implementa a funcionalidade da malha.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
 #include <sys/resource.h>
 #include <sys/time.h>
#include <math.h>

#include "mesh.h"
#include "definitions.h"
#include "helpers.h"



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
    if (! (mesh->x = calloc( (IMAX+1) , sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->y = calloc( (IMAX+1) , sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->xRef = calloc( (IMAX+1) , sizeof(double*)) ) )
    {
        return false;
    }
    if (! (mesh->yRef = calloc( (IMAX+1) , sizeof(double*)) ) )
    {
        return false;
    }
    
    // Alocar memória em y.
    for (int i = 0; i < IMAX+1 ; i++)
    {
        if (! (mesh->x[i] = calloc ( JMAX , sizeof(double))  ) ) 
        {
            return false;
        }
        if (! (mesh->y[i] = calloc ( JMAX , sizeof(double)) ) )
        { 
            return false;
        }
        
        if (! (mesh->xRef[i] = calloc ( JMAX , sizeof(double)) ) )
        { 
            return false;
        }
        if (! (mesh->yRef[i] = calloc ( JMAX , sizeof(double)) ) )
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
    
    for (int i = 0; i < IMAX+1 ; i++)
    {
        free( mesh->x[i] );
        free( mesh->y[i] );
        free( mesh->xRef[i] );
        free( mesh->yRef[i] );
    }
    
    free(mesh->x);
    free(mesh->y);
    free(mesh->xRef);
    free(mesh->yRef);
    return true;
}


bool parabolicMesh(domain * mesh)
{

    for (int j = 1; j < JMAX-1; j++ )
    {
        /** Obtain referenceGrid for J+1 */
        referenceGrid( mesh, j);
        
        /** Solve tridiagonal for x and y */
        
        stepMesh( mesh , j);

        
        
        
        
    }
    
    return true;
}

bool stepMesh ( domain * mesh , int J)
{
    double x_csi, y_csi, x_eta, y_eta;
    
    double psi = 0.0;
    double phi = 0.0;
    
    double * A = malloc ( IMAX * sizeof(double));
    double * B = malloc ( IMAX * sizeof(double));
    double * C = malloc ( IMAX * sizeof(double));
    double * D = malloc ( IMAX * sizeof(double));    
        
    for ( int i = 1; i <= IMAX-2; i++)    
    {
        x_csi = 0.5 * ( mesh->xRef[i+1][J+1] - mesh->xRef[i-1][J+1] );
        y_csi = 0.5 * ( mesh->yRef[i+1][J+1] - mesh->yRef[i-1][J+1] );
        
        x_eta = ( mesh->xRef[i][J+1] - mesh->xRef[i][J] );
        y_eta = ( mesh->yRef[i][J+1] - mesh->yRef[i][J] );
        
        A[i] = x_eta * x_eta + y_eta * y_eta;
        
        B[i] = x_csi * x_eta + y_eta * y_csi;
        
        C[i] = x_csi * x_csi + y_csi * y_csi;
        
        D[i] = (x_csi*y_eta - x_eta*y_csi) * (x_csi*y_eta - x_eta*y_csi);
    }
    
    /** Last and first node */
    x_csi = 0.5 * ( mesh->xRef[1][J+1] - mesh->xRef[IMAX-2][J+1] );
    y_csi = 0.5 * ( mesh->yRef[1][J+1] - mesh->yRef[IMAX-2][J+1] );
    
    x_eta = ( mesh->xRef[IMAX-1][J+1] - mesh->xRef[IMAX-1][J] );
    y_eta = ( mesh->yRef[IMAX-1][J+1] - mesh->yRef[IMAX-1][J] );
    
    A[0] = x_eta * x_eta + y_eta * y_eta;
    A[IMAX-1] = x_eta * x_eta + y_eta * y_eta;
        
    B[0] = x_csi * x_eta + y_eta * y_csi;
    B[IMAX-1] = x_csi * x_eta + y_eta * y_csi;
        
    C[0] = x_csi * x_csi + y_csi * y_csi;
    C[IMAX-1] = x_csi * x_csi + y_csi * y_csi;
        
    D[0] = (x_csi*y_eta - x_eta*y_csi) * (x_csi*y_eta - x_eta*y_csi);
    D[IMAX-1] = (x_csi*y_eta - x_eta*y_csi) * (x_csi*y_eta - x_eta*y_csi);
      
     
    double * lower = malloc( IMAX * sizeof(double));
    double * upper = malloc( IMAX * sizeof(double));
    double * mid = malloc( IMAX * sizeof(double));
    double * d = malloc( IMAX * sizeof(double));   
    
    for (int i = 0; i <= IMAX-1; i++)
    {
        lower[i] = 2.0*A[i] * (1.0 - phi/2.0);
        mid[i] = - 4.0 * (A[i] + C[i]);
        upper[i] = 2.0*A[i] *(1 + phi/2.0);
    }    
     
    for (int i = 1; i <= IMAX-2; i++)
    {    
        d[i] = B[i] * ( mesh->xRef[i+1][J+1] - mesh->x[i+1][J-1] - mesh->xRef[i-1][J+1] + mesh->x[i-1][J-1])
                     - 2.0 * C[i] * (1.0 + psi/2.0) * mesh->xRef[i][J+1]
                     - 2.0 * C[i] * (1.0 - psi/2.0) * mesh->x[i][J-1];
    }
    /** Periodic BC */
    
    d[IMAX-1] = B[IMAX-1] * ( mesh->xRef[1][J+1] - mesh->xRef[1][J+1] - mesh->xRef[IMAX-2][J+1] + mesh->x[IMAX-2][J-1] )
                        - 2.0 * C[IMAX-1] * (1.0 + psi/2.0) * mesh->xRef[IMAX-1][J+1]
                        - 2.0 * C[IMAX-1] * (1.0 - psi/2.0) * mesh->x[IMAX-1][J-1];
    
    d[0] = d[IMAX-1];
    
    solvePeriodicTridiag( lower, mid, upper, d, mesh->x, IMAX , J);
    mesh->x[IMAX-1][J] = 0.5 * (mesh->x[0][J] + mesh->x[IMAX-1][J]);
    mesh->x[0][J] = mesh->x[IMAX-1][J];
    
    
    for (int i = 1; i <= IMAX-2; i++)
    {

    d[i] = B[i] * ( mesh->yRef[i+1][J+1] - mesh->y[i+1][J-1] - mesh->yRef[i-1][J+1] + mesh->y[i-1][J-1])
                 - 2.0 * C[i] * (1 + psi/2.0) * mesh->yRef[i][J+1]
                 - 2.0 * C[i] * (1 - psi/2.0) * mesh->y[i][J-1];
    }
    
    d[IMAX-1] = B[IMAX-1] * ( mesh->yRef[1][J+1] - mesh->y[1][J-1] - mesh->yRef[IMAX-2][J+1] + mesh->y[IMAX-2][J-1])
                 - 2.0 * C[IMAX-1] * (1 + psi/2.0) * mesh->yRef[IMAX-1][J+1]
                 - 2.0 * C[IMAX-1] * (1 - psi/2.0) * mesh->y[IMAX-1][J-1];
    
    /** Periodic BC */
     d[0] = d[IMAX-1];
    
    
  solvePeriodicTridiag( lower, mid, upper, d, mesh->y, IMAX , J);
  
  mesh->y[IMAX-1][J] = 0.5*(mesh->y[0][J]+mesh->y[IMAX-1][J]);
  mesh->y[0][J] = mesh->y[IMAX-1][J];
  
  free(A);
  free(B);
  free(C);
  free(D);
  free(lower);
  free(upper);
  free(mid);
  free(d);        
  
  return true;
  
}

bool innerBC(domain * mesh)
{
    /** Mesh over which configuration? */
    meshCreate( mesh );
    
    return true;
}

bool outerBC(domain * mesh)
{
    /** Create points in a circle with equispaced angles, with radiius R = R_FARFIELD*/
    int NJ = JMAX - 1;
    int NI = IMAX -1;
    double theta;
    
    for (int i = 0; i <= NI; i++)
    {
        theta = ( 2.0*PI / NI ) * i;
        mesh->x[i][NJ] = (double)CHORD/2.0 + R_FARFIELD * (cos(theta)); 
        mesh->y[i][NJ] = - R_FARFIELD * (sin(theta));
    }
    
    return true;
}

bool meshCreate( domain * mesh )
{
    
    int NI = IMAX - 1;
    int np = IMAX;
    double dx; 
    double q;   /** razão da PG */
    double imaxQuarter = ((double) NI/4.0);
    double exponent = NI / 4.0 ;
    q = 1.0 + MESH_STRETCHING;
    
    dx = 0.5 * ( (q - 1.0) / (pow(q,exponent) - 1.0) ) ; /** Força um ponto em 0.5 */
   // dx = (double)CHORD/( (double)(IMAX-1.0)/2.0 );
    double x_centerline[np];
    double y_centerline[np];

    
    /** Generate centerline */
    x_centerline[0] = CHORD;
    y_centerline[0] = 0.0;
    y_centerline[NI] = 0.0;
    
    for (int i = 1; i < IMAX+1; i++)
    {
        
        /** Right bottom */ 
        if (i <= imaxQuarter)
        {
            x_centerline[i] = x_centerline[i-1] - dx ;
            y_centerline[i] = -1.0 * airfoilConfig( x_centerline[i] );
            
            dx = dx*q;

        }
        /** Left bottom */
        if (i > imaxQuarter && i <= 2*imaxQuarter)
        {
            dx = dx/q;
            x_centerline[i] = x_centerline[i-1] - dx;
            y_centerline[i] = -1.0 * airfoilConfig( x_centerline[i] );

        }
        /** Top Left */ 
        if (i > 2*imaxQuarter && i <= 3*imaxQuarter)
        {
            if (i == 2*imaxQuarter+1)
            dx = 0.5 * ( (q - 1.0) / (pow(q,exponent) - 1.0) ) ;        

            x_centerline[i] = x_centerline[i-1] + dx;
            y_centerline[i] = 1.0 * airfoilConfig( x_centerline[i] );
            dx = dx*q;
        }
        /** Top Right */
        if (i > 3*imaxQuarter && i <= 4*imaxQuarter+1)
        {
            dx = dx/q;
            x_centerline[i] = x_centerline[i-1] + dx;
            y_centerline[i] = 1.0 * airfoilConfig( x_centerline[i] );

        }
        
    }
    
    for (int i = 0; i < IMAX; i++)
    {
        mesh->x[i][0] =  x_centerline[i];
        mesh->y[i][0] =  y_centerline[i];
    }
    
    
    /** Periodic Mesh */
    mesh->x[IMAX-1][0] = mesh->x[0][0];
    mesh->y[IMAX-1][0] = mesh->y[0][0];
    
    return true;
}

bool referenceGrid(domain * mesh, int J)
{
    int NJ = JMAX-1;
    int NI = IMAX-1;
    double S_eta, deltaS, R, squareRoot, x_eta, y_eta, x_csi, y_csi, xOrt, yOrt, xInterp, yInterp;
   
    double * s = malloc ( JMAX * sizeof(double) );
        
         /** Calculate Delta S, a function that deals with the streching in the normal direction */
        s[0] = 0;
        s[1] = 1.0 * ( (Q_ESTIRAMENTO - 1.0) / (pow((double)Q_ESTIRAMENTO, NJ) - 1.0) ); /** 1.0 = sum of PG */

        for ( int j = 1; j <= NJ-1; j++)
        {    
            s[j+1] = s[j] + (Q_ESTIRAMENTO) * (s[j] - s[j-1]);
            //s[j+1] = (Q_ESTIRAMENTO) * (s[j]);
            
        }
        
        //deltaS = s[J];
        deltaS = (s[J] - s[J-1]) / (s[NJ] - s[J-1]); 

/** You have to find the reference mesh at J and J+1. */

/** Reference grid at J */

       /** Periodic BC */
        
        /** Run through i */
        for ( int i = 1; i <= IMAX-2; i++) 
        {
            /** Calculate R - distance between j-th line and farfield */
            R = sqrt (  (mesh->x[i][NJ] - mesh->x[i][J-1])*(mesh->x[i][NJ] - mesh->x[i][J-1])
                      + (mesh->y[i][NJ] - mesh->y[i][J-1])*(mesh->y[i][NJ] - mesh->y[i][J-1]));
            /** Calculate S_eta, the derivative of the arc-lenght of a line csi, with respect to eta. */
            S_eta = deltaS * R;
            
            x_csi = 0.5 * ( mesh->x[i+1][J-1] - mesh->x[i-1][J-1]);
            y_csi = 0.5 * (mesh ->y[i+1][J-1] - mesh->y[i-1][J-1]);
            
            /** Calculate x_eta and y_eta s.t. B = 0 (orthogonal) and 1/J increases. */
            squareRoot = sqrt ( x_csi*x_csi + y_csi*y_csi );
                              
            x_eta = - ( y_csi * S_eta ) / ( squareRoot );
            y_eta = + ( x_csi * S_eta ) / ( squareRoot );
      
            /** Calculate Orthogonal Points */
            xOrt = mesh->x[i][J-1] + x_eta;
            yOrt = mesh->y[i][J-1] + y_eta;
            
            /** Calculate Interpolation points */
            xInterp = mesh->x[i][J-1] + deltaS * (mesh->x[i][NJ] - mesh->x[i][J-1]);
            yInterp = mesh->y[i][J-1] + deltaS * (mesh->y[i][NJ] - mesh->y[i][J-1]);
            
            /** Calculate Local Reference Mesh */
            double epsEstiramento =  ((double)J - 1.0) / ((double)NJ - 1.0); 
            
            mesh->xRef[i][J] = epsEstiramento * xInterp + ( 1 - epsEstiramento ) * xOrt;
            mesh->yRef[i][J] = epsEstiramento * yInterp + ( 1 - epsEstiramento ) * yOrt;
        }
    
    /** At Periodic Boundary */

    R = sqrt (  (mesh->x[NI][NJ] - mesh->x[NI][J-1])*(mesh->x[NI][NJ] - mesh->x[NI][J-1])
              + (mesh->y[NI][NJ] - mesh->y[NI][J-1])*(mesh->y[NI][NJ] - mesh->y[NI][J-1]));
    S_eta = deltaS * R;
    x_csi = 0.5 * ( mesh->x[1][J-1] - mesh->x[NI-1][J-1]);
    y_csi = 0.5 * ( mesh->y[1][J-1] - mesh->y[NI-1][J-1]);
    squareRoot = sqrt ( x_csi*x_csi + y_csi*y_csi );
    x_eta = - ( y_csi * S_eta ) / ( squareRoot );
    y_eta = + ( x_csi * S_eta ) / ( squareRoot );
    xOrt = mesh->x[NI][J-1] + x_eta;
    yOrt = mesh->y[NI][J-1] + y_eta;     
    xInterp = mesh->x[NI][J-1] + deltaS * (mesh->x[NI][NJ] - mesh->x[NI][J-1]);
    yInterp = mesh->y[NI][J-1] + deltaS * (mesh->y[NI][NJ] - mesh->y[NI][J-1]);
    double epsEstiramento =  ((double)J - 1.0) / ((double)NJ - 1.0);                     
    mesh->xRef[NI][J] = epsEstiramento * xInterp + ( 1 - epsEstiramento ) * xOrt;
    mesh->yRef[NI][J] = epsEstiramento * yInterp + ( 1 - epsEstiramento ) * yOrt;        

    mesh->xRef[0][J] = mesh->xRef[NI][J];
    mesh->yRef[0][J] = mesh->yRef[NI][J];



/** Reference grid at J+1 */

         /** Periodic BC */
        
        /** Run through i */
        for ( int i = 1; i < IMAX; i++) 
        {
            
            
            /** Calculate R - distance between j-th line and farfield */
            R = sqrt (  (mesh->x[i][NJ] - mesh->xRef[i][J])*(mesh->x[i][NJ] - mesh->xRef[i][J]) 
                      + (mesh->y[i][NJ] - mesh->yRef[i][J])*(mesh->y[i][NJ] - mesh->yRef[i][J]) );
            
            /** Calculate S_eta, the derivative of the arc-lenght of a line csi, with respect to eta. */
            S_eta = deltaS * R;
            
            x_csi = 0.5 * ( mesh->xRef[i+1][J] - mesh->xRef[i-1][J]);
            y_csi = 0.5 * (mesh ->yRef[i+1][J] - mesh->yRef[i-1][J]);
            
            /** Calculate x_eta and y_eta s.t. B = 0 (orthogonal) and 1/J increases. */
            squareRoot = sqrt ( x_csi*x_csi + y_csi*y_csi  );
                              
            x_eta = - ( y_csi * S_eta ) / ( squareRoot );
            y_eta = + ( x_csi * S_eta ) / ( squareRoot );
      
            /** Calculate Orthogonal Points */
            xOrt = mesh->xRef[i][J] + x_eta;
            yOrt = mesh->yRef[i][J] + y_eta;
            
            /** Calculate Interpolation points */
            xInterp = mesh->xRef[i][J] + deltaS * (mesh->x[i][NJ] - mesh->xRef[i][J]);
            yInterp = mesh->yRef[i][J] + deltaS * (mesh->y[i][NJ] - mesh->yRef[i][J]);
            
            /** Calculate Local Reference Mesh */
            double epsEstiramento =  ((double)(J+1) - 1.0) / ((double)NJ - 1.0); 
            
            mesh->xRef[i][J+1] = epsEstiramento * xInterp + ( 1 - epsEstiramento ) * xOrt;
            mesh->yRef[i][J+1] = epsEstiramento * yInterp + ( 1 - epsEstiramento ) * yOrt;
        }
    
    /** At Periodic Boundary */
    
    R = sqrt (  (mesh->x[NI][NJ] - mesh->xRef[NI][J])*(mesh->x[NI][NJ] - mesh->xRef[NI][J]) 
              + (mesh->y[NI][NJ] - mesh->yRef[NI][J])*(mesh->y[NI][NJ] - mesh->yRef[NI][J]) );
    S_eta = deltaS * R;
    x_csi = 0.5 * ( mesh->xRef[1][J] - mesh->xRef[NI-1][J]);
    y_csi = 0.5 * ( mesh->yRef[1][J] - mesh->yRef[NI-1][J]);
    squareRoot = sqrt ( x_csi*x_csi + y_csi*y_csi  );
    x_eta = - ( y_csi * S_eta ) / ( squareRoot );
    y_eta = + ( x_csi * S_eta ) / ( squareRoot );            
    xOrt = mesh->xRef[NI][J] + x_eta;
    yOrt = mesh->yRef[NI][J] + y_eta;    
    xInterp = mesh->xRef[NI][J] + deltaS * (mesh->x[NI][NJ] - mesh->xRef[NI][J]);
    yInterp = mesh->yRef[NI][J] + deltaS * (mesh->y[NI][NJ] - mesh->yRef[NI][J]);
    epsEstiramento =  ((double)(J+1) - 1.0) / ((double)NJ - 1.0); 
            
    mesh->xRef[NI][J+1] = epsEstiramento * xInterp + ( 1 - epsEstiramento ) * xOrt;
    mesh->yRef[NI][J+1] = epsEstiramento * yInterp + ( 1 - epsEstiramento ) * yOrt;    
    
    mesh->xRef[0][J+1] = mesh->xRef[NI][J+1];
    mesh->yRef[0][J+1] = mesh->yRef[NI][J+1];
    
    free(s);
    
    return true;
}

double airfoilConfig(double x)
{
    return 5.0*th*CHORD*(0.2969*sqrt(fabs(x)) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1015*x*x*x*x);
    
    return 2.0 * th * x * ( 1.0 - x ) ;
}

