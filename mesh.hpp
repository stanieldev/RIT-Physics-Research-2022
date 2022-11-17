/*
 * File:	mesh.hpp
 * Author:	Stanley Goodwin
 * Stores the simulation's mesh characteristics.
 */
#pragma once
#ifndef MESH_H
#define MESH_H

#include "droplet.hpp"
#include "substrate.hpp"
#include "node.hpp"

#define MAX_POINTS 122


struct Mesh {
    Droplet current_droplet;
    Substrate current_surface;
    double volume = 0;
    double pressure = 0;
    double gamma = 0;  // Gamma factor
private:
    int res;   // Resolution of the mesh (MUST BE ODD)
    int res1;  // Resolution - 1
    int res2;  // Half the resolution floored to nearest integer
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
    Node(*swap_nodes)[MAX_POINTS];

public:
    Mesh();
    Mesh(int _resolution, Droplet _droplet, Substrate _surface);

    void change_mesh_resolution(int _new_size);
    void initialize(double _initial_contact_angle_degrees);
    void iterate(int _iteration_count);
    double calc_volume();
    double calc_pressure();

    // File print functions of mesh
    void fprint_nodes();     // Current nodes to file
    void fprint_volume();    // Current volume to file
    void fprint_pressure();  // Current pressure to file
    void fprint_gamma();     // Current gamma factor to file

    // Console print functions of mesh
    void cprint_nodes();     // Current nodes to console
    void cprint_volume();    // Current volume to console
    void cprint_pressure();  // Current pressure to console
    void cprint_gamma();     // Current gamma factor to console

private:
    void _swap_nodes();

    void _initialize(double _initial_contact_angle);
    void _iterate(int _iteration_count);
    void _update_volume();
    void _update_pressure();

    double _contact_angle(int i, int j);
    double _approximate_surface(Node n1, Node n2, Node n3, Node n4, Node n5, Node n6, Node selected);

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