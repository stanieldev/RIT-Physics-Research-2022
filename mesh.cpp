/*
 * File:    mesh.cpp
 * Author:  Stanley Goodwin
 * Mesh-related functions.
 */
#include <iostream>
#include <fstream>
#include <string>
#include "fmath.hpp"
#include "timing.h"
#include "mesh.h"



/***********************************************************
**                  Simulation Functions                  **
***********************************************************/



/*
 * Initializes the mesh node array to prepare it for iteration.
 * Creates an initial mesh in the shape of a spherical cap.
 * https://en.wikipedia.org/wiki/Spherical_coordinate_system
 * @brief   Initialize mesh's nodes with a spherical cap.
 */
void Mesh::initialize(double initial_contact_angle)
{
    // Constants
    const double k = 1.0 / sin(initial_contact_angle);  // Radius scale factor
    const double h = -k * cos(initial_contact_angle);   // Vertical shift (spherical cap)

    // Spherical coordinates
    double θ, Δθ;  // Polar angle
    double φ, Δφ;  // Azimuthal angle

    // Derived coordinates (spherical -> cartesian where ρ = 1)
    double z, r;  // The azimuth [cosine, sine] component of the droplet contact radius
    double x, y;  // The  polar  [cosine, sine] component of the azimuth radius

    // Other misc
    double volume_factor;

    // Initialization start
    std::cout << "Calculating Initial Node Mesh... ";
    time start = hrc::now();

    // Defining the azimuthal angle
    Δφ = initial_contact_angle / m_res2;
    φ = 0;

    // Create middle point
    m_current_nodes[m_res2][m_res2] = Node(0, 0, k * cos(φ) + h);
    φ += Δφ;

    // Loop over radial segments of the spherical cap surface
    for (int j = 1; j <= m_res2; j++)  // Radial incrementer
    {
        // Defining variables
        r = k * sin(φ);
        z = k * cos(φ) + h;

        // Defining the polar angle
        Δθ = (PI / 2) / (2 * j);
        θ = 0;

        // Loop over polar rotations about the z-axis
        for (int i = 0; i < 2 * j; i++)  // Angular incrementer
        {
            // Defining variables
            x = r * cos(θ);
            y = r * sin(θ);

            // Initialize node position vectors
            m_current_nodes[m_res2 - j + i][m_res2 - j] = Node( x,  y, z);  // Bottom
            m_current_nodes[m_res2 + j][m_res2 - j + i] = Node(-y,  x, z);  // Right
            m_current_nodes[m_res2 + j - i][m_res2 + j] = Node(-x, -y, z);  // Top
            m_current_nodes[m_res2 - j][m_res2 + j - i] = Node( y, -x, z);  // Left

            // Increment polar angle
            θ += Δθ;
        }

        // Increment azimuthal angle
        φ += Δφ;
    }

    // Scale node heights for a mesh volume ~ 1.1 * expected volume
    volume_factor = 1.1 * droplet.vkappa / volume();
    for (int i = 0; i < m_res; i++)  // Current angle index
    {
        for (int j = 0; j < m_res; j++)  // Current radius index
        {
            m_current_nodes[i][j].z *= 1.1 * volume_factor;
        }
    }

    // Calculate mesh characteristics
    _volume = volume();
    _pressure = pressure();

    // Print step
    current_iter = 0; 
    iterations_run = 0;  
    fprint_nodes();
    iterations_run = 1;

    // Initialization conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << _duration_string(start, stop) << ")\n";
}

/*
 * Iterates the mesh iteration_count steps toward final geometry.
 * @brief   Iterate mesh's node arrays.
 */
void Mesh::iterate(int iteration_count)
{
    // Constants
    double k1 = surface.kp;
    double k2 = surface.kg;

    // Variables
    Node mean;  // The mean value of nodes
    Node diff;  // The change in the nodes
    double _contact_angle;  // The nodes's contact angle
    current_iter = 1;


    // Iteration start
    std::cout << "Iterating Mesh... ";
    time start = hrc::now();

    // Iterate the mesh toward the expected volume
    for (int λ = 1; λ <= iteration_count; λ++)
    {
        // Swap current & previous node memory addresses
        m_swap_nodes = m_previous_nodes;
        m_previous_nodes = m_current_nodes;
        m_current_nodes = m_swap_nodes;

        // Create current nodes from previous nodes
        for (int i = 0; i < m_res; i++)  // Current angle index
        {
            for (int j = 0; j < m_res; j++)  // Current radius index
            {
                // Boundary points
                if (i == 0 || i == m_res1 || j == 0 || j == m_res1) {

                    // Find the node's contact angle
                    _contact_angle = contact_angle(i, j);

                    // If node is on the printed region(s) and under the printed slip angle
                    if (surface.slips_on_printed(m_previous_nodes[i][j], _contact_angle))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i][j + 1]) / 3;
                            diff = m_previous_nodes[i + 1][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j + 1].z * k1);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i + 1][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j - 1].z * k1);
                        }

                        // Top-right boundary point
                        else if (i == m_res1 && j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i - 1][j] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i - 1][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j - 1].z * k1);
                        }

                        // Bottom-right boundary point
                        else if (i == m_res1 && j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i - 1][j] + m_previous_nodes[i][j + 1]) / 3;
                            diff = m_previous_nodes[i - 1][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j + 1].z * k1);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i - 1][j]) / 3;
                            diff = m_previous_nodes[i][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i][j + 1].z * k1);
                        }

                        // Right side
                        else if (i == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i][j + 1] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i - 1][j] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j].z * k1);
                        }

                        // Top side
                        else if (j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i - 1][j]) / 3;
                            diff = m_previous_nodes[i][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i][j - 1].z * k1);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i][j + 1] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i + 1][j] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j].z * k1);
                        }

                        // Create new node
                        m_current_nodes[i][j] = mean * α + diff * β;
                        m_current_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If node is on the printed region(s) and under the gap region slip angle
                    else if (surface.slips_on_surface(m_previous_nodes[i][j], _contact_angle))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i][j + 1]) / 3;
                            diff = m_previous_nodes[i + 1][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j + 1].z * k2);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i + 1][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j - 1].z * k2);
                        }

                        // Top-right boundary point
                        else if (i == m_res1 && j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i - 1][j] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i - 1][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j - 1].z * k2);
                        }

                        // Bottom-right boundary point
                        else if (i == m_res1 && j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i - 1][j] + m_previous_nodes[i][j + 1]) / 3;
                            diff = m_previous_nodes[i - 1][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j + 1].z * k2);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i - 1][j]) / 3;
                            diff = m_previous_nodes[i][j + 1] + vector_gradient(i, j) * (m_previous_nodes[i][j + 1].z * k2);
                        }

                        // Right side
                        else if (i == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i][j + 1] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i - 1][j] + vector_gradient(i, j) * (m_previous_nodes[i - 1][j].z * k2);
                        }

                        // Top side
                        else if (j == m_res1)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i + 1][j] + m_previous_nodes[i - 1][j]) / 3;
                            diff = m_previous_nodes[i][j - 1] + vector_gradient(i, j) * (m_previous_nodes[i][j - 1].z * k2);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (m_previous_nodes[i][j] + m_previous_nodes[i][j + 1] + m_previous_nodes[i][j - 1]) / 3;
                            diff = m_previous_nodes[i + 1][j] + vector_gradient(i, j) * (m_previous_nodes[i + 1][j].z * k2);
                        }

                        // Create new node
                        m_current_nodes[i][j] = mean * α + diff * β;
                        m_current_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If the node's contact angle is steeper than the slip angles (pinned)
                    else {
                        m_current_nodes[i][j] = m_previous_nodes[i][j];
                    }
                }

                // Internal points
                else
                {
                    m_current_nodes[i][j] = m_previous_nodes[i][j] + NetNormalForce(i, j) + NetTangentialForce(i, j);
                }
            }
        }

        // Print step
        if (λ % node_print_interval == 0) { fprint_nodes(); }
        current_iter++;

        // Calculate mesh characteristics
        _volume = volume();
        _pressure = pressure();
    }

    // Print debug info if enabled
    current_iter = 0;
    fprint_nodes();
    iterations_run++;

    // Iteration conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << _duration_string(start, stop) << ")\n";
}



/***********************************************************
**                Characteristic Functions                **
***********************************************************/

/*
 * Returns the approximate current volume of the mesh surface.
 * @brief	Current mesh volume.
 * @return	volume	double	The volume of the current surface.
 */
double Mesh::volume()
{
    // Variables
    double V = 0.0;   // Volume accumulator
    double dA = 0.0;  // Small area segment
    double h = 0.0;   // Mean height of dA
    Node v1, v2, v3, v4;

    // Sum all the volume segments
    for (int i = 0; i < m_res1; i++)
    {
        for (int j = 0; j < m_res1; j++)
        {
            // Assign variables
            v1 = m_current_nodes[i + 1][  j  ] - m_current_nodes[  i  ][  j  ];
            v2 = m_current_nodes[  i  ][j + 1] - m_current_nodes[  i  ][  j  ];
            v3 = m_current_nodes[i + 1][j + 1] - m_current_nodes[  i  ][j + 1];
            v4 = m_current_nodes[i + 1][j + 1] - m_current_nodes[i + 1][  j  ];

            // Calculate characteristics
            h = (m_current_nodes[i][j] + m_current_nodes[i + 1][j] + m_current_nodes[i][j + 1] + m_current_nodes[i + 1][j + 1]).z;
            dA = (v1.x) * (v2.y) - (v2.x) * (v1.y) + (v3.x) * (v4.y) - (v4.x) * (v3.y);
            
            // Add new volume to total volume
            V += h * dA;
        }
    }

    // Return volume
    return V / 8.0;
}

/*
 * Returns the approximate current pressure of the mesh surface.
 * @brief	Current mesh pressure.
 * @return	presure	double	The pressure of the current surface.
 */
double Mesh::pressure()
{
    // Variables
    double _volume = volume();
    double _pressure;
    double base;
    double expo;

    // Iterate Gamma factor
    _gamma += δ * droplet.radius3 * (droplet.vkappa - _volume);

    // Calculate pressure
    base = droplet.vkappa / (_volume * _volume);
    expo = (1. + 0.1 * (droplet.vkappa / _volume + _volume / droplet.vkappa - 2)) / 3.0;
    _pressure = exp(_gamma) * (σ / µ) * pow(base, expo);

    // Return pressure
    return _pressure;
}







// Status: INCOMPLETE
// TODO: Same as VectorGradient.
/*
 * Returns the surface contact angle at point (i,j).
 * Note: Node(i, j) is only valid for boundary points.
 *
 * @brief	The contact angle of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The contact angle IN RADIANS.
 */
double Mesh::contact_angle(int i, int j)
{
    // Variables
    Node p1, p2, p3, p4, p5, p6;
    const double step = 0.1 / m_res;

    // Bottom-left edge point
    if (i == 0 && j == 0)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i + 1][j];
        p3 = m_previous_nodes[i][j + 1];
        p4 = m_previous_nodes[i + 1][j + 1];
        p5 = m_previous_nodes[i + 2][j + 1];
        p6 = m_previous_nodes[i + 1][j + 2];
    }

    // Top-left edge point
    else if (i == 0 && j == m_res1)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i + 1][j];
        p3 = m_previous_nodes[i][j - 1];
        p4 = m_previous_nodes[i + 1][j - 1];
        p5 = m_previous_nodes[i + 2][j - 1];
        p6 = m_previous_nodes[i + 1][j - 2];
    }

    // Top-right edge point
    else if (i == m_res1 && j == m_res1)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i - 1][j];
        p3 = m_previous_nodes[i][j - 1];
        p4 = m_previous_nodes[i - 1][j - 1];
        p5 = m_previous_nodes[i - 2][j - 1];
        p6 = m_previous_nodes[i - 1][j - 2];
    }

    // Bottom-right edge point
    else if (i == m_res1 && j == 0)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i - 1][j];
        p3 = m_previous_nodes[i][j + 1];
        p4 = m_previous_nodes[i - 1][j + 1];
        p5 = m_previous_nodes[i - 2][j + 1];
        p6 = m_previous_nodes[i - 1][j + 2];
    }

    // Bottom edge
    else if (j == 0)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i - 1][j];
        p3 = m_previous_nodes[i + 1][j];
        p4 = m_previous_nodes[i][j + 1];
        p5 = m_previous_nodes[i - 1][j + 1];
        p6 = m_previous_nodes[i + 1][j + 1];
    }

    // Right edge
    else if (i == m_res1)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i][j - 1];
        p3 = m_previous_nodes[i][j + 1];
        p4 = m_previous_nodes[i - 1][j];
        p5 = m_previous_nodes[i - 1][j - 1];
        p6 = m_previous_nodes[i - 1][j + 1];
    }

    // Top edge
    else if (j == m_res1)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i - 1][j];
        p3 = m_previous_nodes[i + 1][j];
        p4 = m_previous_nodes[i][j - 1];
        p5 = m_previous_nodes[i - 1][j - 1];
        p6 = m_previous_nodes[i + 1][j - 1];
    }

    // Left edge
    else if (i == 0)
    {
        p1 = m_previous_nodes[i][j];
        p2 = m_previous_nodes[i][j - 1];
        p3 = m_previous_nodes[i][j + 1];
        p4 = m_previous_nodes[i + 1][j];
        p5 = m_previous_nodes[i + 1][j - 1];
        p6 = m_previous_nodes[i + 1][j + 1];
    }


    // Generate an approximating surface using those 6 nodes
    double matrix[6][6] = {
        { p1.x * p1.x, p1.x * p1.y, p1.y * p1.y, p1.x, p1.y, 1.0 },
        { p2.x * p2.x, p2.x * p2.y, p2.y * p2.y, p2.x, p2.y, 1.0 },
        { p3.x * p3.x, p3.x * p3.y, p3.y * p3.y, p3.x, p3.y, 1.0 },
        { p4.x * p4.x, p4.x * p4.y, p4.y * p4.y, p4.x, p4.y, 1.0 },
        { p5.x * p5.x, p5.x * p5.y, p5.y * p5.y, p5.x, p5.y, 1.0 },
        { p6.x * p6.x, p6.x * p6.y, p6.y * p6.y, p6.x, p6.y, 1.0 },
    };
    double output[6] = { p1.z, p2.z, p3.z, p4.z, p5.z, p6.z };
    double coeff[6];

    // Determine the coefficients of the polynomial
    cramer(matrix, output, coeff);


    // Take a step in the gradient direction
    int sign = (p1.x * p1.x + p1.y * p1.y < p4.x* p4.x + p4.y * p4.y) ? 1 : -1;
    Node new_node = p1 - vector_gradient(i, j) * step * sign;

    // Find the new node's z-component using the polynomial approximation
    double new_z = coeff[0] * new_node.x * new_node.x +
        coeff[1] * new_node.x * new_node.y +
        coeff[2] * new_node.y * new_node.y +
        coeff[3] * new_node.x +
        coeff[4] * new_node.y +
        coeff[5];

    // Return the inverse tangent of the slope of the change in the node position
    return atan2(new_z, -step * sign);
}



/***********************************************************
**                    Vector Functions                    **
***********************************************************/

// Status: INCOMPLETE
// TODO: Direction of gradient not taken into account (of node above it is closer or further from origin).
/*
 * Returns the normalized gradient vector at point (i,j).
 * Note: Node(i, j) is only valid for boundary points.
 *
 * @brief	The normalized gradient vector of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The normalized gradient.
 */
Node Mesh::vector_gradient(int i, int j)
{
    // Variables
    Node v;
    int n = m_res1;
    int sign = 1;

    // Choose vector to rotate by 90 degrees
         if (i == 0 && j == 0) { v = m_previous_nodes[i + 1][j] - m_previous_nodes[i][j + 1]; }  // Bottom-left edge point
    else if (i == 0 && j == n) { v = m_previous_nodes[i][j - 1] - m_previous_nodes[i + 1][j]; }  // Top-left edge point
    else if (i == n && j == n) { v = m_previous_nodes[i - 1][j] - m_previous_nodes[i][j - 1]; }  // Top-right edge point
    else if (i == n && j == 0) { v = m_previous_nodes[i][j + 1] - m_previous_nodes[i - 1][j]; }  // Bottom-right edge point
    else if (j == 0)           { v = m_previous_nodes[i + 1][j] - m_previous_nodes[i - 1][j]; }  // Bottom edge
    else if (i == n)           { v = m_previous_nodes[i][j + 1] - m_previous_nodes[i][j - 1]; }  // Right edge
    else if (j == n)           { v = m_previous_nodes[i - 1][j] - m_previous_nodes[i + 1][j]; }  // Top edge
    else if (i == 0)           { v = m_previous_nodes[i][j - 1] - m_previous_nodes[i][j + 1]; }  // Left edge

    // Return normalized gradient vector
    return Node(v.y, -v.x, 0).normalize();
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the normal vector at point(i, j) of the surface.
 * Note: Node(i, j) is only valid for non-boundary points.
 *
 * @brief	Normal vector at surface.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The normal vector at that point.
 */
Node Mesh::vector_normal(int i, int j)
{
    // Variables
    Node v1, v2;

    /*// Top-left diagonal
    if (i == m_res1 - j && i < m_res2)
    {
        v1 = m_previous_nodes[i - 1][j + 1] - m_previous_nodes[i + 1][j - 1];  // Radius change
        v2 = m_previous_nodes[i][j - 1] - m_previous_nodes[i + 1][j];          // Polar angle change
    }

    // Bottom-right diagonal
    else if (i == m_res1 - j && i > m_res2)
    {
        v1 = m_previous_nodes[i + 1][j - 1] - m_previous_nodes[i - 1][j + 1];  // Radius change
        v2 = m_previous_nodes[i][j + 1] - m_previous_nodes[i - 1][j];          // Polar angle change
    }

    // Bottom-left diagonal
    else if (i == j && i < m_res2)
    {
        v1 = m_previous_nodes[i - 1][j - 1] - m_previous_nodes[i + 1][j + 1];  // Radius change
        v2 = m_previous_nodes[i + 1][j] - m_previous_nodes[i][j + 1];          // Polar angle change
    }

    // Top-right diagonal
    else if (i == j && i > m_res2)
    {
        v1 = m_previous_nodes[i + 1][j + 1] - m_previous_nodes[i - 1][j - 1];  // Radius change
        v2 = m_previous_nodes[i - 1][j] - m_previous_nodes[i][j - 1];          // Polar angle change
    }

    // Non-diagonal points & center point
    else {
        v1 = m_previous_nodes[i + 1][j] - m_previous_nodes[i - 1][j];  // Normal neighbor differential
        v2 = m_previous_nodes[i][j + 1] - m_previous_nodes[i][j - 1];  // Normal neighbor differential
    }*/

    // Generate intersecting vectors
    v1 = m_previous_nodes[i + 1][j] - m_previous_nodes[i - 1][j];
    v2 = m_previous_nodes[i][j + 1] - m_previous_nodes[i][j - 1];

    // Return normal vector
    return cross_product(v1, v2).normalize();
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the tangential vector to the surface centered at point (i,j).
 *
 * @brief	Mesh tangential vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The tangent vector at that point.
 */
Node Mesh::vector_tangent(int i, int j)
{
    // Variables
    Node v_cross = m_previous_nodes[i][j - 1] + m_previous_nodes[i][j + 1] + m_previous_nodes[i + 1][j] + m_previous_nodes[i - 1][j] - m_previous_nodes[i][j] * 4;
    Node v_normal = vector_normal(i, j);

    // Return tangent vector
    return v_cross - v_cross.proj(v_normal);
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the mean curvature at point (i,j).
 *
 * @brief	Mean curvature vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The mean curvature vector at that point.
 */
Node Mesh::vector_mean_curvature(int i, int j)
{
    // Variables
    Node v_i = m_previous_nodes[i + 1][j] - m_previous_nodes[i - 1][j];  // 2nd Difference in X
    Node v_j = m_previous_nodes[i][j + 1] - m_previous_nodes[i][j - 1];  // 2nd Difference in Y

    Node v_b = (m_previous_nodes[i][j - 1] - m_previous_nodes[i][j]).normalize();  // Bottom
    Node v_l = (m_previous_nodes[i - 1][j] - m_previous_nodes[i][j]).normalize();  // Left
    Node v_r = (m_previous_nodes[i + 1][j] - m_previous_nodes[i][j]).normalize();  // Right
    Node v_t = (m_previous_nodes[i][j + 1] - m_previous_nodes[i][j]).normalize();  // Left

    Node s1 = (v_l + v_r) * v_j.det();
    Node s2 = (v_t + v_b) * v_i.det();

    // Return mean curvature
    return (s1 + s2) * 0.5;  // Second derivative approximation
}



/***********************************************************
**                 Force Vector Functions                 **
***********************************************************/

// Status: COMPLETE & VERIFIED
/**
 * Returns the mesh curvature force at point (i,j).
 *
 * @brief	Mesh curvature force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The curvature force at that point.
 */
Node Mesh::CurvatureForce(int i, int j)
{
    // Variables
    Node normal_vector = vector_normal(i, j);

    // Return mean curvature normal vector
    return vector_mean_curvature(i, j).proj(normal_vector) * (σ / µ);
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the pressure force at point (i,j).
 *
 * @brief	Mesh pressure force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The pressure force vector at that point.
 */
Node Mesh::PressureForce(int i, int j)
{
    // Variables
    Node v1 = m_previous_nodes[i - 1][j] - m_previous_nodes[i][j - 1];
    Node v2 = m_previous_nodes[i - 1][j] - m_previous_nodes[i][j + 1];
    Node v3 = m_previous_nodes[i + 1][j] - m_previous_nodes[i][j - 1];
    Node v4 = m_previous_nodes[i + 1][j] - m_previous_nodes[i][j + 1];

    // Calculated cross products
    Node vec1 = cross_product(v1, v2);
    Node vec2 = cross_product(v3, v4);

    // Return normal pressure force
    return vector_normal(i, j) * (m_pressure * 0.25) * (vec1.det() + vec2.det());
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the tangential force to the surface centered at point (i,j).
 *
 * @brief	Mesh tangential force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The tangent force at that point.
 */
Node Mesh::TangentialForce(int i, int j)
{
    return vector_tangent(i, j) * (τ / µ);
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the net normal force to the surface centered at point (i,j).
 *
 * @brief	Mesh net normal force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The net normal force at that point.
 */
Node Mesh::NetNormalForce(int i, int j)
{
    return PressureForce(i, j) + CurvatureForce(i, j);
}

// Status: COMPLETE & VERIFIED
/**
 * Returns the net tangential force to the surface centered at point (i,j).
 *
 * @brief	Mesh net tangential force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The net tangent force at that point.
 */
Node Mesh::NetTangentialForce(int i, int j)
{
    return TangentialForce(i, j);
}
