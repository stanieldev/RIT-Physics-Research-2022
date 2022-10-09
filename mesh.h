//
// File: mesh.h
// Author: Stanley Goodwin
// Creation Date: 9/17/2021
// Last Modified: 6/16/2022
// Credit to Kara Maki for skeleton code.
//
#pragma once
#ifndef MESH_H
#define MESH_H

#include <string>
#include "node.h"
#include "droplet.h"

#define MAX_POINTS 51
#define MAX_TIMES 2500


/*
The class that defines the fishnet mesh that will take the shape of the droplet.
*/
class Mesh {
private:
    // Control panel (for now)
    int _total_iterations = 2100;   // Total mesh iterations
    int _resolution = 41;   // The count of angular & radial subdivisions
    int _res1 = _resolution - 1;

    // Mesh constants
    const double Δt = 5.E-8;   // delta time unit
    const double τ = 1 * 1.E2; // Tangential constant (tau)
    const double µ = 4.e-5;    // Frictional constant (mu)
    const double σ = 60;       // Coefficient of surface tention (sigma)
    const double α = 0.5;      // boundary relaxation factor (coding alpha)
    const double β = 1 - α;    // boundary relaxation factor (coding beta)
    const double δ = 1.0e4;    // relaxation factor for the pressure

    // Arrays
    Node* _node_array[MAX_POINTS][MAX_POINTS];
    Node* _normal_vector[MAX_POINTS][MAX_POINTS];



    // Figure this shit out
protected:
    double m_dP = 1.0E5;
    double m_PressureFactor = 0;

    double _volume[MAX_TIMES];
    double _pressure[MAX_TIMES];




public:
    void CInitNodes();
    void CBoundary(int λ);
    void CInterior(int λ);





    // Normal

    // Vector
    void Iterate();
    double vCurrentVolume(int count);
    double vPressure(int count);

    Node NormalVectorBottom(int i, int count);
    Node NormalVectorRight(int j, int count);
    Node NormalVectorTop(int i, int count);
    Node NormalVectorLeft(int j, int count);






    

    Node NetNormalForce(int i, int j, int count);
    Node NetTangentialForce(int i, int j, int count);
    

    
    double PressureForce(int i, int j, int count);
    Node MeanCurvature(int i, int j, int count);
    double SigmaForce(int i, int j, int count);
    Node NormalVector(int i, int j, int count);
    
    Node TangentPart(int i, int j, int count);





    

    void vprintCurrentMassInformation(std::string s, Droplet droplet);

    Mesh();

};


#endif MESH_H