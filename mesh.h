//
// File: mesh.h
// Author: Stanley Goodwin
// Creation Date: 9/17/2021
// Last Modified: 7/4/2022
// Credit to Kara Maki for skeleton code
//
#pragma once
#ifndef MESH_H
#define MESH_H

#include <vector>
#include "droplet.h"
#include "node.h"

#define PI 3.141593265358979323846
#define MAX_POINTS 122


/*
The class that defines the fishnet mesh that will take the shape of the droplet.
*/
class Mesh {
private:
    // User parameters
    int _nλ;    // Total mesh iterations
    int _res;   // The number of angular & radial subdivisions
    int _res1;  // ^^^ minus 1 for simpification

    // Strength parameters
    double µ = 8.0e2;    // Frictional decay constant
    double σ = 60.0;     // Coefficient of surface tension
    double τ = 1.0e2;    // Tangential "spring" constant
    double δ = 1.0e4;    // relaxation factor for the pressure
    double α = 0.5;      // Boundary relaxation factor (coding alpha)
    double β = 1.0 - α;  // Conjugate relaxation factor (coding beta)

    // Initial droplet constants
    const double θ = PI / 4;           // Initial contact angle
    const double θ_c = 75 * PI / 180;  // Printed region contact angle

    // Simple variables
    Droplet _droplet;
    double _volume = NULL;    // Mesh's current volume
    double _pressure = NULL;  // Mesh's current pressure
    double _Γ = 0;

    // Array variables
    Node _nodes_prev[MAX_POINTS][MAX_POINTS];
    Node _nodes_curr[MAX_POINTS][MAX_POINTS];

    // Pointer variables
    Node(*_prev_nodes)[MAX_POINTS] = _nodes_prev;
    Node(*_curr_nodes)[MAX_POINTS] = _nodes_curr;
    Node(*_swap_nodes)[MAX_POINTS] = NULL;


    // Functions
    bool OnPrintedRegion(int i, int j);  // Tests if node is on printed region

    double Volume();    // Returns the mesh's current volume
    double Pressure();  // Returns the mesh's current pressure

    Node NormalVector(int i, int j);

    Node TangentPart(int i, int j);
    Node TangentialForce(int i, int j);
    Node NetTangentialForce(int i, int j);

    Node MeanCurvatureIntegral(int i, int j);
    Node CurvatureForce(int i, int j);
    Node PressureForce(int i, int j);
    Node NetNormalForce(int i, int j);


public:
    // Public functions
    void InitializeNodes();    // Create the initial mesh "guess"
    void Iterate();            // Iterate the mesh to final droplet shape
    void PrintCurrent(int λ);  // Prints the current iteration's nodes to files

    // Constructor & Destructor
    Mesh(Droplet droplet, int resolution, int total_iteration_count);
};


#endif MESH_H