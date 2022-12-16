/*
 * File:	mesh.hpp
 * Author:	Stanley Goodwin
 * Stores the simulation's mesh characteristics.
 */
#pragma once
#ifndef MESH_H
#define MESH_H

#include <string>
#include "droplet.hpp"
#include "substrate.hpp"
#include "node.hpp"

#define MAX_POINTS 122


struct Mesh {
    Droplet current_droplet;
    Substrate current_surface;
    double volume = 0.0;
    double pressure = 0.0;
    double gamma = 0.0;
    bool print_nodes_bool = false;
    std::string output_folder;

private:
    int res;   // Resolution of the mesh (MUST BE ODD)
    int res1;  // Resolution - 1
    int res2;  // Half the resolution floored to nearest integer
    int iterations_complete = 0;
    double grad_step;
    double µ = 8.0e2;    // Frictional decay constant (8.0e2)
    double σ = 6.0e1;    // Coefficient of surface tension
    double τ = 1.0e2;    // Tangential "spring" constant
    double δ = 1.0e7;    // relaxation factor for the pressure
    double α = 0.5;      // Boundary relaxation factor (coding alpha)
    double β = 1.0 - α;  // Conjugate relaxation factor (coding beta)
    Node nodes_a[MAX_POINTS][MAX_POINTS];
    Node nodes_b[MAX_POINTS][MAX_POINTS];
    Node(*current_nodes)[MAX_POINTS] = nodes_a;
    Node(*previous_nodes)[MAX_POINTS] = nodes_b;
    Node(*swap_nodes)[MAX_POINTS] = NULL;

public:
    Mesh();
    Mesh(int _resolution, Droplet _droplet, Substrate _surface);
    void initialize(double _initial_contact_angle_degrees);
    void iterate(int _iteration_count);
    void print_nodes();

private:
    void _initialize(double _initial_contact_angle);
    void _iterate(int _iteration_count);
    double _calculate_volume();
    double _calculate_pressure();
    double _update_volume();
    double _update_pressure();

    double _contact_angle(int i, int j);

    Node _vector_gradient(int i, int j);
    Node _vector_normal(int i, int j);
    Node _vector_mean_curvature(int i, int j);
    Node _vector_tangent(int i, int j);

    Node _force_curvature(int i, int j);
    Node _force_pressure(int i, int j);
    Node _force_tangential(int i, int j);
    Node _force_net_normal(int i, int j);
    Node _force_net_tangential(int i, int j);
};

#endif MESH_H