/**
 * @file	mesh.cpp
 * @brief	The mesh function definitions.
 *
 * @author	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *          sfg5318@g.rit.edu
 *
 * Creation Date: 6/16/2022
 * Last Modified: 7/25/2022
 */
#include <iostream>
#include <fstream>
#include <string>
#include "fmath.h"
#include "misc.h"
#include "mesh.h"
#include "new_substrate.h"



/***********************************************************
**                  Simulation Functions                  **
***********************************************************/

// Status: INCOMPLETE
// TODO: Allocate nodes to heap instead of memory
Mesh::Mesh()
{
    // Heap allocation
    /*Node** _prev_array;
    _prev_array = new Node*[_res];
    for (int i = 0; i < _res; i++)
    {
        _prev_array[i] = new Node[_res];
    }*/
}

/**
 * Function description:
 *   Initializes the mesh node array to prepare it for iteration.
 *   Creates an initial mesh in the shape of a spherical cap.
 * 
 * Resources:
 *   https://en.wikipedia.org/wiki/Spherical_coordinate_system
 */
// Status: COMPLETE & VERIFIED
void Mesh::initialize(double θi)
{
    // Constants
    const double k = 1 / sin(θi);   // Radius scale factor
    const double h = -k * cos(θi);  // Vertical shift (spherical cap)

    // Spherical coordinates
    double θ, Δθ;  // Polar angle
    double φ, Δφ;  // Azimuthal angle

    // Derived coordinates (spherical -> cartesian where ρ = 1)
    double z, r;  // The azimuth [cosine, sine] component of the droplet contact radius
    double x, y;  // The  polar  [cosine, sine] component of the azimuth radius


    // Initialization start
    std::cout << "Calculating Initial Node Mesh... ";
    time start = hrc::now();
    
    // Defining the azimuthal angle
    Δφ = θi / m_res2;
    φ = 0;

    // Create middle point
    m_curr_nodes[m_res2][m_res2] = Node(0, 0, k * cos(φ) + h);
    φ += Δφ;
    
    // Loop over radial segments of the spherical cap surface
    for (int j = 1; j <= m_res2; j++)  // Radial incrementer
    {
        // Defining variables
        r = k * sin(φ);
        z = (k * cos(φ) + h);

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
            m_curr_nodes[m_res2 - j + i][m_res2 - j    ] = Node( x,  y, z);  // Bottom
            m_curr_nodes[m_res2 + j    ][m_res2 - j + i] = Node(-y,  x, z);  // Right
            m_curr_nodes[m_res2 + j - i][m_res2 + j    ] = Node(-x, -y, z);  // Top
            m_curr_nodes[m_res2 - j    ][m_res2 + j - i] = Node( y, -x, z);  // Left

            // Increment polar angle
            θ += Δθ;
        }

        // Increment azimuthal angle
        φ += Δφ;
    }

    // Scale node heights for a mesh volume ~ 1.1 * expected volume
    const double s = 1.1 * drop_κ / volume();
    for (int i = 0; i < m_res; i++)  // Current angle index
    {
        for (int j = 0; j < m_res; j++)  // Current radius index
        {
            m_curr_nodes[i][j].z *= s;
        }
    }

    // Calculate mesh characteristics
    m_iteration = 0;
    m_volume = volume();
    m_pressure = pressure();

    // Initialization conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << _duration_string(start, stop) << ")\n";


    // Save current nodes to file
    fprint_nodes();
}

/**
 * Function description:
 *   Iterates the mesh iteration_count steps toward final geometry.
 */
// Status: INCOMPLETE
// TODO: Change For Loop into While Loop (using volumes)
void Mesh::iterate(int iteration_count, Substrate surface)
{
    // Constants
    const double k1 = surface.k1;
    const double k2 = surface.k2;

    // Variables
    Node mean;  // The mean value of nodes
    Node diff;  // The change in the nodes
    double _contact_angle;
    bool _on_printed_region;


    // Iteration start
    std::cout << "Iterating Mesh... ";
    time start = hrc::now();

    // Iterate the mesh toward the expected volume
    for (int λ = 1; λ <= iteration_count; λ++)
    {
        // Swap current & previous node memory addresses
        m_swap_nodes = m_prev_nodes;
        m_prev_nodes = m_curr_nodes;
        m_curr_nodes = m_swap_nodes;

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
                    if (surface.slips_on_printed(m_prev_nodes[i][j], _contact_angle))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i][j + 1]) / 3;
                            diff = m_prev_nodes[i + 1][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j + 1].z * k1);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i + 1][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j - 1].z * k1);
                        }

                        // Top-right boundary point
                        else if (i == m_res1 && j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i - 1][j] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i - 1][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j - 1].z * k1);
                        }

                        // Bottom-right boundary point
                        else if (i == m_res1 && j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i - 1][j] + m_prev_nodes[i][j + 1]) / 3;
                            diff = m_prev_nodes[i - 1][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j + 1].z * k1);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i - 1][j]) / 3;
                            diff = m_prev_nodes[i][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i][j + 1].z * k1);
                        }

                        // Right side
                        else if (i == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i][j + 1] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i - 1][j] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j].z * k1);
                        }

                        // Top side
                        else if (j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i - 1][j]) / 3;
                            diff = m_prev_nodes[i][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i][j - 1].z * k1);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i][j + 1] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i + 1][j] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j].z * k1);
                        }

                        // Create new node
                        m_curr_nodes[i][j] = mean * α + diff * β;
                        m_curr_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If node is on the printed region(s) and under the gap region slip angle
                    /*else if (_contact_angle < θ_d && !_on_printed_region)
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i][j + 1]) / 3;
                            diff = m_prev_nodes[i + 1][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j + 1].z * k2);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i + 1][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j - 1].z * k2);
                        }

                        // Top-right boundary point
                        else if (i == m_res1 && j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i - 1][j] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i - 1][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j - 1].z * k2);
                        }

                        // Bottom-right boundary point
                        else if (i == m_res1 && j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i - 1][j] + m_prev_nodes[i][j + 1]) / 3;
                            diff = m_prev_nodes[i - 1][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j + 1].z * k2);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i - 1][j]) / 3;
                            diff = m_prev_nodes[i][j + 1] + vector_gradient(i, j) * (m_prev_nodes[i][j + 1].z * k2);
                        }

                        // Right side
                        else if (i == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i][j + 1] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i - 1][j] + vector_gradient(i, j) * (m_prev_nodes[i - 1][j].z * k2);
                        }

                        // Top side
                        else if (j == m_res1)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i + 1][j] + m_prev_nodes[i - 1][j]) / 3;
                            diff = m_prev_nodes[i][j - 1] + vector_gradient(i, j) * (m_prev_nodes[i][j - 1].z * k2);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (m_prev_nodes[i][j] + m_prev_nodes[i][j + 1] + m_prev_nodes[i][j - 1]) / 3;
                            diff = m_prev_nodes[i + 1][j] + vector_gradient(i, j) * (m_prev_nodes[i + 1][j].z * k2);
                        }

                        // Create new node
                        m_curr_nodes[i][j] = mean * α + diff * β;
                        m_curr_nodes[i][j].z = 0;  // Boundary condition
                    }*/

                    // If the node's contact angle is steeper than the slip angles (pinned)
                    else {
                        m_curr_nodes[i][j] = m_prev_nodes[i][j];
                    }
                }

                // Internal points
                else
                {
                    m_curr_nodes[i][j] = m_prev_nodes[i][j] + NetNormalForce(i, j) + NetTangentialForce(i, j);
                }
            }
        }

        // Calculate mesh characteristics
        m_volume = volume();
        m_pressure = pressure();
    }

    // Add 1 to volume iteration counter
    m_iteration++;

    // Iteration conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << _duration_string(start, stop) << ")\n";


    // Save current nodes to file
    fprint_nodes();
}



/***********************************************************
**                Characteristic Functions                **
***********************************************************/

/**
 * Returns the approximate current volume of the mesh surface.
 *
 * @brief	Current mesh volume.
 * @return	volume	double	The volume of the current surface.
 */
// Status: COMPLETE & VERIFIED
double Mesh::volume()
{
    // Variables
    double V = 0;   // Volume accumulator
    double dA = 0;  // Small area segment
    double h = 0;   // Mean height of dA
    Node v1, v2, v3, v4;

    // Sum all the volume segments
    for (int i = 0; i < m_res1; i++)
    {
        for (int j = 0; j < m_res1; j++)
        {
            // Assign variables
            v1 = m_curr_nodes[i + 1][  j  ] - m_curr_nodes[  i  ][  j  ];
            v2 = m_curr_nodes[  i  ][j + 1] - m_curr_nodes[  i  ][  j  ];
            v3 = m_curr_nodes[i + 1][j + 1] - m_curr_nodes[  i  ][j + 1];
            v4 = m_curr_nodes[i + 1][j + 1] - m_curr_nodes[i + 1][  j  ];

            // Calculate characteristics
            h = (m_curr_nodes[i][j] + m_curr_nodes[i + 1][j] + m_curr_nodes[i][j + 1] + m_curr_nodes[i + 1][j + 1]).z;
            dA = (v1.x) * (v2.y) - (v2.x) * (v1.y) + (v3.x) * (v4.y) - (v4.x) * (v3.y);
            
            // Add new volume to total volume
            V += h * dA;
        }
    }

    // Return volume
    return V / 8;
}

/**
 * Returns the approximate current pressure of the mesh surface.
 *
 * @brief	Current mesh pressure.
 * @return	presure	double	The pressure of the current surface.
 */
// Status: INCOMPLETE (90%)
// TODO: Calls volume twice and Find better convergence for gamma.
double Mesh::pressure()
{
    // Variables
    double _volume = volume();  // Current Volume
    double _pressure;
    double base;
    double expo;

    // Iterate Gamma factor
    m_gamma += δ * drop_r3 * (drop_κ - _volume);

    // Calculate pressure
    base = drop_κ / (_volume * _volume);
    expo = (1. + 0.1 * (drop_κ / _volume + _volume / drop_κ - 2)) / 3.0;
    _pressure = exp(m_gamma) * (σ / µ) * pow(base, expo);

    // Return pressure
    return _pressure;
}



/***********************************************************
**                   Boundary Functions                   **
***********************************************************/

/*
 * Returns the normalized gradient vector at point (i,j).
 * Note: Node(i, j) is only valid for boundary points.
 *
 * @brief	The normalized gradient vector of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The normalized gradient.
 */
// Status: INCOMPLETE
// TODO: Direction of gradient not taken into account (of node above it is closer or further from origin).
Node Mesh::vector_gradient(int i, int j)
{
    // Variables
    Node v;
    int sign = 1;

    // Bottom-left edge point
    if (i == 0 && j == 0) 
    {
        v = m_prev_nodes[i + 1][j] - m_prev_nodes[i][j + 1];
    }

    // Top-left edge point
    else if (i == 0 && j == m_res1) 
    { 
        v = m_prev_nodes[i][j - 1] - m_prev_nodes[i + 1][j];
    }

    // Top-right edge point
    else if (i == m_res1 && j == m_res1) 
    { 
        v = m_prev_nodes[i - 1][j] - m_prev_nodes[i][j - 1]; 
    }

    // Bottom-right edge point
    else if (i == m_res1 && j == 0) 
    { 
        v = m_prev_nodes[i][j + 1] - m_prev_nodes[i - 1][j]; 
    }

    // Bottom edge
    else if (j == 0) 
    { 
        v = m_prev_nodes[i + 1][j] - m_prev_nodes[i - 1][j]; 
    }

    // Right edge
    else if (i == m_res1) 
    { 
        v = m_prev_nodes[i][j + 1] - m_prev_nodes[i][j - 1]; 
    }

    // Top edge
    else if (j == m_res1) 
    { 
        v = m_prev_nodes[i - 1][j] - m_prev_nodes[i + 1][j]; 
    }

    // Left edge
    else if (i == 0) 
    { 
        v = m_prev_nodes[i][j - 1] - m_prev_nodes[i][j + 1]; 
    }

    // Return normalized gradient vector
    return Node(v.y, -v.x, 0).normalize();
}

/*
 * Returns the surface contact angle at point (i,j).
 * Note: Node(i, j) is only valid for boundary points.
 *
 * @brief	The contact angle of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The contact angle IN RADIANS.
 */
// Status: INCOMPLETE
// TODO: Same as VectorGradient.
double Mesh::contact_angle(int i, int j)
{
    // Variables
    Node p1, p2, p3, p4, p5, p6;
    const double step = 0.1 / m_res;

    // Bottom-left edge point
    if (i == 0 && j == 0) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i + 1][j];
        p3 = m_prev_nodes[i][j + 1];
        p4 = m_prev_nodes[i + 1][j + 1];
        p5 = m_prev_nodes[i + 2][j + 1];
        p6 = m_prev_nodes[i + 1][j + 2];
    }

    // Top-left edge point
    else if (i == 0 && j == m_res1) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i + 1][j];
        p3 = m_prev_nodes[i][j - 1];
        p4 = m_prev_nodes[i + 1][j - 1];
        p5 = m_prev_nodes[i + 2][j - 1];
        p6 = m_prev_nodes[i + 1][j - 2];
    }

    // Top-right edge point
    else if (i == m_res1 && j == m_res1) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i - 1][j];
        p3 = m_prev_nodes[i][j - 1];
        p4 = m_prev_nodes[i - 1][j - 1];
        p5 = m_prev_nodes[i - 2][j - 1];
        p6 = m_prev_nodes[i - 1][j - 2];
    }

    // Bottom-right edge point
    else if (i == m_res1 && j == 0) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i - 1][j];
        p3 = m_prev_nodes[i][j + 1];
        p4 = m_prev_nodes[i - 1][j + 1];
        p5 = m_prev_nodes[i - 2][j + 1];
        p6 = m_prev_nodes[i - 1][j + 2];
    }

    // Bottom edge
    else if (j == 0) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i - 1][j];
        p3 = m_prev_nodes[i + 1][j];
        p4 = m_prev_nodes[i][j + 1];
        p5 = m_prev_nodes[i - 1][j + 1];
        p6 = m_prev_nodes[i + 1][j + 1];
    }

    // Right edge
    else if (i == m_res1) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i][j - 1];
        p3 = m_prev_nodes[i][j + 1];
        p4 = m_prev_nodes[i - 1][j];
        p5 = m_prev_nodes[i - 1][j - 1];
        p6 = m_prev_nodes[i - 1][j + 1];
    }

    // Top edge
    else if (j == m_res1) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i - 1][j];
        p3 = m_prev_nodes[i + 1][j];
        p4 = m_prev_nodes[i][j - 1];
        p5 = m_prev_nodes[i - 1][j - 1];
        p6 = m_prev_nodes[i + 1][j - 1];
    }

    // Left edge
    else if (i == 0) 
    {
        p1 = m_prev_nodes[i][j];
        p2 = m_prev_nodes[i][j - 1];
        p3 = m_prev_nodes[i][j + 1];
        p4 = m_prev_nodes[i + 1][j];
        p5 = m_prev_nodes[i + 1][j - 1];
        p6 = m_prev_nodes[i + 1][j + 1];
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
    double output[6] = {p1.z, p2.z, p3.z, p4.z, p5.z, p6.z};
    double coeff[6];

    // Determine the coefficients of the polynomial
    cramer(matrix, output, coeff);


    
    // Take a step in the gradient direction
    int sign = (p1.x * p1.x + p1.y * p1.y < p4.x * p4.x + p4.y * p4.y) ? 1 : -1;
    Node new_node = p1 - vector_gradient(i, j) * step * sign;
    
    // Find the new node's z-component using the polynomial approximation
    double new_z = coeff[0] * new_node.x * new_node.x +
                   coeff[1] * new_node.x * new_node.y +
                   coeff[2] * new_node.y * new_node.y +
                   coeff[3] * new_node.x +
                   coeff[4] * new_node.y +
                   coeff[5];
    
    // Return the inverse tangent of the slope of the change in the node position
    return atan2(new_z, - step * sign);
}

/*
 * Returns whether a point (i,j) is on any printed region.
 *
 * @brief	If the node(i,j) is on printed region.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	onreg	bool	True if the node is on the printed region.
 *                          False if the node is not on the printed region.
 */
// Status: INCOMPLETE
// TODO: Make more user friendly.
//bool Mesh::on_printed_region(int i, int j)
//{
//    // Absolute positions
//    double abs_y_val = abs(m_prev_nodes[i][j].y) - 0.5 * w_p;
//
//    // Region booleans
//    bool r1 = (    -0.5 * w_p    <= abs_y_val && abs_y_val <= 0 * w_g + 0 * w_p);
//    bool r2 = (1 * w_g + 0 * w_p <= abs_y_val && abs_y_val <= 1 * w_g + 1 * w_p);
//    /*bool r3 = (2 * w_g + 1 * w_p <= abs_y_val && abs_y_val <= 2 * w_g + 2 * w_p);
//    bool r4 = (3 * w_g + 2 * w_p <= abs_y_val && abs_y_val <= 3 * w_g + 3 * w_p);
//    bool r5 = (4 * w_g + 3 * w_p <= abs_y_val && abs_y_val <= 4 * w_g + 4 * w_p);
//    bool r6 = (5 * w_g + 4 * w_p <= abs_y_val && abs_y_val <= 5 * w_g + 5 * w_p);
//    bool r7 = (6 * w_g + 5 * w_p <= abs_y_val && abs_y_val <= 6 * w_g + 6 * w_p);
//    bool r8 = (7 * w_g + 6 * w_p <= abs_y_val && abs_y_val <= 7 * w_g + 7 * w_p);*/
//
//    // Return if point is on printed region
//    return r1 || r2; //  || r3 || r4 || r5 || r6 || r7 || r8
//}



/***********************************************************
**                   Interior Functions                   **
***********************************************************/

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

    // Top-left diagonal
    if (i == m_res1 - j && i < m_res2)
    {  
        v1 = m_prev_nodes[i - 1][j + 1] - m_prev_nodes[i + 1][j - 1];  // Radius change
        v2 = m_prev_nodes[i][j - 1] - m_prev_nodes[i + 1][j];          // Polar angle change
    }

    // Bottom-right diagonal
    else if (i == m_res1 - j && i > m_res2)
    {  
        v1 = m_prev_nodes[i + 1][j - 1] - m_prev_nodes[i - 1][j + 1];  // Radius change
        v2 = m_prev_nodes[i][j + 1] - m_prev_nodes[i - 1][j];          // Polar angle change
    }

    // Bottom-left diagonal
    else if (i == j && i < m_res2)
    {  
        v1 = m_prev_nodes[i - 1][j - 1] - m_prev_nodes[i + 1][j + 1];  // Radius change
        v2 = m_prev_nodes[i + 1][j] - m_prev_nodes[i][j + 1];          // Polar angle change
    }

    // Top-right diagonal
    else if (i == j && i > m_res2)
    {
        v1 = m_prev_nodes[i + 1][j + 1] - m_prev_nodes[i - 1][j - 1];  // Radius change
        v2 = m_prev_nodes[i - 1][j] - m_prev_nodes[i][j - 1];          // Polar angle change
    }

    // Non-diagonal points & center point
    else {
        v1 = m_prev_nodes[i + 1][j] - m_prev_nodes[i - 1][j];  // Normal neighbor differential
        v2 = m_prev_nodes[i][j + 1] - m_prev_nodes[i][j - 1];  // Normal neighbor differential
    }

    // Return normal vector
    return cross_product(v1, v2).normalize();
}


/**
 * Returns the mean curvature at point (i,j).
 *
 * @brief	Mean curvature vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The mean curvature vector at that point.
 */
Node Mesh::MeanCurvatureIntegral(int i, int j)
{
    // Variables
    Node v_i = m_prev_nodes[i + 1][j] - m_prev_nodes[i - 1][j];  // 2nd Difference in X
    Node v_j = m_prev_nodes[i][j + 1] - m_prev_nodes[i][j - 1];  // 2nd Difference in Y

    Node v_b = (m_prev_nodes[i][j - 1] - m_prev_nodes[i][j]).normalize();  // Bottom
    Node v_l = (m_prev_nodes[i - 1][j] - m_prev_nodes[i][j]).normalize();  // Left
    Node v_r = (m_prev_nodes[i + 1][j] - m_prev_nodes[i][j]).normalize();  // Right
    Node v_t = (m_prev_nodes[i][j + 1] - m_prev_nodes[i][j]).normalize();  // Left

    Node s1 = (v_l + v_r) * v_j.det();
    Node s2 = (v_t + v_b) * v_i.det();
    
    // Return mean curvature
    return (s1 + s2) / 2;  // Second derivative approximation
}

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
    return MeanCurvatureIntegral(i, j).proj(normal_vector) * (σ / µ);
}

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
    Node v1 = m_prev_nodes[i - 1][j] - m_prev_nodes[i][j - 1];
    Node v2 = m_prev_nodes[i - 1][j] - m_prev_nodes[i][j + 1];
    Node v3 = m_prev_nodes[i + 1][j] - m_prev_nodes[i][j - 1];
    Node v4 = m_prev_nodes[i + 1][j] - m_prev_nodes[i][j + 1];

    // Calculated cross products
    Node vec1 = cross_product(v1, v2);
    Node vec2 = cross_product(v3, v4);

    // Return normal pressure force
    return vector_normal(i, j) * (m_pressure / 4) * (vec1.det() + vec2.det());
}

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


/**
 * Returns the tangential vector to the surface centered at point (i,j).
 *
 * @brief	Mesh tangential vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The tangent vector at that point.
 */
Node Mesh::TangentPart(int i, int j)
{
    // Variables
    Node v_cross = m_prev_nodes[i][j - 1] + m_prev_nodes[i][j + 1] + m_prev_nodes[i + 1][j] + m_prev_nodes[i - 1][j] - m_prev_nodes[i][j] * 4;
    Node v_normal = vector_normal(i, j);

    // Return vector
    return v_cross - v_cross.proj(v_normal);
}

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
    return TangentPart(i, j) * (τ / µ);
}

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
