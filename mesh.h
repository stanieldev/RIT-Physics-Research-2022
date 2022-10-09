//
// File: mesh.h
// Author: Stanley Goodwin
// Creation Date: 9/17/2021
// Last Modified: 6/23/2022
// Credit to Kara Maki for skeleton code.
//
#pragma once
#ifndef MESH_H
#define MESH_H

#include <string>
#include "node.h"
#include "droplet.h"

#define MAX_POINTS 101
#define MAX_TIMES 2500


/*
The class that defines the fishnet mesh that will take the shape of the droplet.
*/
class Mesh {
private:
    // Control panel (for now)
    const int _total_iterations = 1200;  // Total mesh iterations
    const int _resolution = 21;          // The count of angular & radial subdivisions

    // Mesh strength parameters
    const double µ = 8.0e2;  // Frictional decay constant
    const double σ = 60.0;   // Coefficient of surface tension
    const double τ = 1.0e2;  // Tangential "spring" constant
    const double δ = 1.0e4;  // relaxation factor for the pressure
    const double α = 0.5;    // Boundary relaxation factor (coding alpha)
    const double β = 1 - α;  // Conjugate relaxation factor (coding beta)

    // Mesh constants
    const int _res1 = _resolution - 1;
    
    // Mesh variables
    Node* _node_array[MAX_POINTS][MAX_POINTS];

    Node _prev_nodes;
    Node _curr_nodes;

    double _volume;
    double _pressure;
    double _Γ = 0;
    

    // General functions
    bool OnPrintedRegion(int i, int j, int λ);  // Finds if node is on printed region
    
    // Normal Vector functions
    Node NormalVectorBottom(int i, int count);
    Node NormalVectorRight(int j, int count);
    Node NormalVectorTop(int i, int count);
    Node NormalVectorLeft(int j, int count);
    Node NormalVector(int i, int j, int count);
    
    // Calculation functions
    double Volume(int λ);    // Returns the mesh's volume
    double Pressure(int λ);  // Returns the mesh's pressure

    Node TangentPart(int i, int j, int count);
    Node TangentialForce(int i, int j, int count);
    Node NetTangentialForce(int i, int j, int count);

    Node MeanCurvatureIntegral(int i, int j, int count);
    Node CurvatureForce(int i, int j, int count);
    Node PressureForce(int i, int j, int count);
    Node NetNormalForce(int i, int j, int count);
    

public:
    // Public Functions
    void InitializeNodes();  // Create the initial mesh "guess"
    void Iterate();          // Iterate the mesh to final droplet shape
    void vprintCurrentMassInformation(std::string s, Droplet droplet);
    Mesh();
};


#endif MESH_H