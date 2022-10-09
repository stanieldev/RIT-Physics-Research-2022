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
#ifndef MESH_H
#define MESH_H

#include <vector>
#include "node.h"
#include "settings.h"
#include "new_substrate.h"

#define MAX_POINTS 122


/*
 * Mesh characteristics class.
 *
 * Stores all the necessary data of the fishnet mesh,
 * as well as the initialization and simulation of 
 * iterations on the mesh in order to match best the
 * droplet in question.
 */
class Mesh: MeshSettings {
public:

    // Simulation constants
    const int m_res  = 121;        // The number of radial subdivisions (MUST BE ODD)
    const int m_res1 = m_res - 1;  // One less than m_res
    const int m_res2 = m_res / 2;  // Floored m_res / 2

    // Droplet settings
    const double drop_radius = 1.9407025E-3;  // The contact radius of the droplet
    const double drop_volume = 3.0E-9;        // Expected volume of the droplet
    const double drop_r3 = drop_radius * drop_radius * drop_radius;
    double drop_κ = drop_volume / drop_r3;

    // Printing variables
    int m_iteration = 0;  // Current volume mesh iteration number

    // Mesh variables
    double m_volume = NULL;    // Current volume
    double m_pressure = NULL;  // Current pressure
    double m_gamma = 0;        // Current gamma factor (m_gamma)

    // Mesh array variables
    Node nodes_array_previous[MAX_POINTS][MAX_POINTS];
    Node nodes_array_current[MAX_POINTS][MAX_POINTS];

    // Mesh pointer variables
    Node(*m_prev_nodes)[MAX_POINTS] = nodes_array_previous;
    Node(*m_curr_nodes)[MAX_POINTS] = nodes_array_current;
    Node(*m_swap_nodes)[MAX_POINTS] = NULL;

public:
    // COMPLETE
    void initialize(double θi);  // Create the initial mesh
    



    // Simulation functions



    
    Mesh();             // Constructor
    void iterate(int iteration_count, Substrate surface);  // Iterate iteration_count steps
    double volume();    // Returns the mesh's current volume

private:
    // COMPLETE
    
    double pressure();  // Returns the mesh's current pressure





    // Characteristic functions
    

    // Boundary functions
    Node vector_gradient(int i, int j);
    double contact_angle(int i, int j);
    //bool on_printed_region(int i, int j);  // Tests if node is on printed region

    // Interior functions
    Node vector_normal(int i, int j);

    Node MeanCurvatureIntegral(int i, int j);
    Node CurvatureForce(int i, int j);
    Node PressureForce(int i, int j);
    Node NetNormalForce(int i, int j);

    Node TangentPart(int i, int j);
    Node TangentialForce(int i, int j);
    Node NetTangentialForce(int i, int j);
    
public:

    // File print functions
    void fprint_nodes();     // Current nodes to file
    void fprint_volume();    // Current volume to file
    void fprint_pressure();  // Current pressure to file
    void fprint_gamma();     // Current gamma factor to file

    // Console print functions
    void cprint_nodes();     // Current nodes to console
    void cprint_volume();    // Current volume to console
    void cprint_pressure();  // Current pressure to console
    void cprint_gamma();     // Current gamma factor to console
};


#endif MESH_H