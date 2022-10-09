//
//  WorkingDrop.h
//  DropletShape
//
//  Created by Kara Maki on 9/17/21.
//  Copyright © 2021 Kara Maki. All rights reserved.
//

#ifndef Redo_h
#define Redo_h

#include <stdio.h>
#include <string>

#define PI 3.141593
#define MAX_POINTS 100
#define MAX_TIMES 10000

class Drop {

protected:
    int m_nPoints;  // Number of points on a curve
    int m_mPoints;    // n - horizontal, m - vertical
    int n1, n2, n3;
    int m1, m2, m3;

    double* m_pdXt[MAX_POINTS][MAX_POINTS];    // we calc them
    double* m_pdYt[MAX_POINTS][MAX_POINTS];
    double* m_pdZt[MAX_POINTS][MAX_POINTS];

    double pdPress[MAX_TIMES];    // an array of pressure
    double pdVol[MAX_TIMES];    // and array of volume
    double theta[MAX_POINTS][MAX_POINTS];

    int m_nTimes;    // number of time points
    int m_nCount;
    double m_dDeltaT;  // delta time unit
    double m_dTan;     // Tangential constant
    double m_dFric;    // Frictional constant
    double m_dMass;    // mass of a point
    double m_dSigma;   // Coefficient of surface tention
    double m_dAlpha;    // boundary relaxation factor
    double m_PressRelax; //relaxation factor for the pressure

    double m_dXmin, m_dXmax,    // phys limits by X
        m_dYmin, m_dYmax;        // phys limits by Y
    double m_dV;                // volume of a drop
    double m_dTheta;            // the angle between a surface and a drop
    double m_dP;
    double m_PressureFactor;

public:
    void CalculateTimeArrays();
    double CurrentVolume(int count);
    double KnownNormalForceX(int j, int k, int count);
    double KnownNormalForceY(int j, int k, int count);
    double KnownNormalForceZ(int j, int k, int count);
    double KnownTangentialForceX(int j, int k, int count);
    double KnownTangentialForceY(int j, int k, int count);
    double KnownTangentialForceZ(int j, int k, int count);
    double NormalVectorBottomX(int i, int count);
    double NormalVectorBottomY(int i, int count);
    double NormalVectorBottomZ(int i, int count);
    double NormalVectorRightX(int j, int count);
    double NormalVectorRightY(int j, int count);
    double NormalVectorRightZ(int j, int count);
    double NormalVectorTopX(int i, int count);
    double NormalVectorTopY(int i, int count);
    double NormalVectorTopZ(int i, int count);
    double NormalVectorLeftX(int j, int count);
    double NormalVectorLeftY(int j, int count);
    double NormalVectorLeftZ(int j, int count);
    double NormalVectorX(int j, int k, int count);
    double NormalVectorY(int j, int k, int count);
    double NormalVectorZ(int j, int k, int count);
    double Pressure(int count);
    double PressureConst(int j, int k, int count);
    double PressureForceX(int j, int k, int count);
    double PressureForceY(int j, int k, int count);
    double PressureForceZ(int j, int k, int count);
    double SigmaDiffX(int j, int k, int count);
    double SigmaDiffY(int j, int k, int count);
    double SigmaDiffZ(int j, int k, int count);
    double SigmaForceX(int j, int k, int count);
    double SigmaForceY(int j, int k, int count);
    double SigmaForceZ(int j, int k, int count);
    double TangentPartX(int j, int k, int count);
    double TangentPartY(int j, int k, int count);
    double TangentPartZ(int j, int k, int count);
    double Theta(int i, int j, int count);
    void Initialization();
    void printCurrentMassInformation(std::string s);

    Drop(); //constructor of the drop class

};

#endif /* Redo_hpp */
