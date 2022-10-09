/**
 * @file	mesh.h
 * @brief	A class that stores mesh variables and functions.
 *
 * @author	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 9/17/2021
 * Last Modified: 7/6/2022
 */
#pragma once
#ifndef NODE_H
#define NODE_H

#include <vector>
#include "node.h"
#include "settings.h"

#define MAX_POINTS 122


/**
 * Mesh characteristics class.
 *
 * Stores all the necessary data of the fishnet mesh,
 * as well as the initialization and simulation of 
 * iterations on the mesh in order to match best the
 * droplet in question.
 */
class Mesh: MeshSettings {
protected:
    // Precalculated
    int _res1 = _res - 1;  // 1 less than resolution
    double half_p = 0.5 * w_p;  // Half the printed region size

    // Simple variables
    double _volume = NULL;    // Mesh's current volume
    double _pressure = NULL;  // Mesh's current pressure
    double _Γ = 0;            // Mesh's gamma factor (pressure)

    // Array variables
    Node _nodes_prev[MAX_POINTS][MAX_POINTS];
    Node _nodes_curr[MAX_POINTS][MAX_POINTS];

    // Pointer variables
    Node(*_prev_nodes)[MAX_POINTS] = _nodes_prev;
    Node(*_curr_nodes)[MAX_POINTS] = _nodes_curr;
    Node(*_swap_nodes)[MAX_POINTS] = NULL;

private:
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

    
    Mesh();  // Constructor
};


#endif MESH_H