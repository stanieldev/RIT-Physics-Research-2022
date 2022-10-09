/**
 * File:	mesh.h
 * Brief:	A class that stores mesh variables and functions.
 *
 * Author:	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 9/17/2021
 * Last Modified: 7/26/2022
 */
#pragma once
#ifndef MESH_H
#define MESH_H

#include "substrate.h"
#include "droplet.h"
#include "node.h"

#define PI 3.141593265358979323846
#define MAX_POINTS 122


/*
 * Mesh characteristics class.
 *
 * Stores all the necessary data of the fishnet mesh,
 * as well as the initialization and simulation of 
 * iterations on the mesh in order to match best the
 * droplet in question.
 */
class Mesh {
public:
    int iterations_run = 0;
    int current_iter = 0;
    double m_volume = 0;    // Current volume
    double m_pressure = 0;  // Current pressure
    double m_gamma = 0;     // Current gamma factor (m_gamma)

private:
    Droplet droplet;
    Substrate surface;

    int node_print_interval;  // How often to save intermediate steps

    int m_res;  // The number of radial subdivisions (MUST BE ODD)
    int m_res1;  // One less than m_res
    int m_res2;  // Floored m_res / 2

    double µ = 8.0e2;    // Frictional decay constant (8.0e2)
    double σ = 6.0e1;    // Coefficient of surface tension
    double τ = 1.0e2;    // Tangential "spring" constant
    double δ = 1.0e7;    // relaxation factor for the pressure
    double α = 0.5;      // Boundary relaxation factor (coding alpha)
    double β = 1.0 - α;  // Conjugate relaxation factor (coding beta)

    Node m_previous_node_array[MAX_POINTS][MAX_POINTS];
    Node m_current_node_array[MAX_POINTS][MAX_POINTS];

    Node(*m_previous_nodes)[MAX_POINTS] = m_previous_node_array;
    Node(*m_current_nodes)[MAX_POINTS] = m_current_node_array;
    Node(*m_swap_nodes)[MAX_POINTS] = NULL;


public:
    Mesh(int _resolution, Droplet _droplet, Substrate _surface, int _node_print_interval)
    {
        m_res = _resolution;
        m_res1 = m_res - 1;
        m_res2 = m_res / 2;
        droplet = _droplet;
        surface = _surface;
        node_print_interval = _node_print_interval;
    };
    
    void rescale(double volume_factor);
    void initialize(double initial_contact_angle);
    void iterate(int iteration_count);

private:
    double volume();    // Returns the mesh's current volume
    double pressure();  // Returns the mesh's current pressure
    double contact_angle(int i, int j);

    Node vector_gradient(int i, int j);
    Node vector_normal(int i, int j);
    Node vector_mean_curvature(int i, int j);
    Node vector_tangent(int i, int j);

    Node CurvatureForce(int i, int j);
    Node PressureForce(int i, int j);
    Node TangentialForce(int i, int j);
    Node NetNormalForce(int i, int j);
    Node NetTangentialForce(int i, int j);
    
    
public:

    //// File print functions
    void fprint_nodes();     // Current nodes to file
    //void fprint_volume();    // Current volume to file
    //void fprint_pressure();  // Current pressure to file
    //void fprint_gamma();     // Current gamma factor to file

    //// Console print functions
    //void cprint_nodes();     // Current nodes to console
    //void cprint_volume();    // Current volume to console
    //void cprint_pressure();  // Current pressure to console
    //void cprint_gamma();     // Current gamma factor to console
};


#endif MESH_H