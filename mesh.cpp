/*
 * File:    mesh.cpp
 * Author:  Stanley Goodwin
 * Stores the simulation's mesh characteristics.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include "timing.hpp"
#include "fmath.hpp"
#include "mesh.hpp"



/*
 * Definition of constructing a mesh object.
 * @brief	Default mesh contructor.
 */
Mesh::Mesh()
{
    res = 25;
    res1 = res - 1;
    res2 = res / 2;
    grad_step = 0.1 / res;
}

/*
 * Definition of constructing a mesh object.
 * @brief	Value-defined mesh contructor.
 * @param	_mesh_resolution    double
 * @param	_droplet            Droplet
 * @param	_surface            Substrate
 */
Mesh::Mesh(int _mesh_resolution, Droplet _droplet, Substrate _surface)
{
    res = _mesh_resolution;
    res1 = _mesh_resolution - 1;
    res2 = _mesh_resolution / 2;
    grad_step = 0.1 / res;
    current_droplet = _droplet;
    current_surface = _surface;
    assert(_mesh_resolution % 2 == 1);  // Resolution must be odd number
};



/***********************************************************
**                     Utility Functions                  **
***********************************************************/

/*
 * Prints the contents of the node arrays to a text file.
 * @brief   Save mesh's nodes to file.
 */
void Mesh::print_nodes()
{
    // Variables
    const std::string file_name = std::to_string(res) + "x" + std::to_string(res);
    const std::string file_iter = "_" + std::to_string(iterations_complete);
    const std::string file_end = ".txt";
    std::ofstream data, file_x, file_y, file_z;
    Node _node;

    // Open folders
    data.open(output_folder + "\\" + file_name + "_data" + file_end, std::ios_base::app);
    file_x.open(output_folder + "\\" + file_name + file_iter + "_x" + file_end);
    file_y.open(output_folder + "\\" + file_name + file_iter + "_y" + file_end);
    file_z.open(output_folder + "\\" + file_name + file_iter + "_z" + file_end);

    // Push data to folders
    data << res << "x" << res << " " << iterations_complete << " " << volume << " " << pressure << " " << gamma << "\n";
    for (int i = 0; i < res; i++)
    {
        for (int j = 0; j < res; j++)
        {
            // Set current node
            _node = current_nodes[i][j] * current_droplet.contact_radius;
        
            // Send data to the files
            file_x << _node.x << " "; file_y << _node.y << " "; file_z << _node.z << " ";
        }
        // Add new line
        file_x << "\n"; file_y << "\n"; file_z << "\n";
    }
        
    // Close file streams
    data.close(); file_x.close(); file_y.close(); file_z.close();

    // Add 1 to iterator
    iterations_complete += 1;
};



/***********************************************************
**                  Simulation Functions                  **
***********************************************************/

/*
 * Initializes the mesh node array to prepare it for iteration.
 * Creates an initial mesh in the shape of a spherical cap.
 * https://en.wikipedia.org/wiki/Spherical_coordinate_system
 * @brief   Initialize mesh's nodes with a spherical cap.
 */
void Mesh::initialize(double _initial_contact_angle_degrees)  // Public
{
    // Initialization start
    std::cout << "Calculating Initial Node Mesh... ";
    time start = hrc::now();

    // Run function
    constexpr auto DEG_TO_RAD = 3.141593265358979323846 / 180.0;
    _initialize(_initial_contact_angle_degrees * DEG_TO_RAD);

    // Print if enabled
    iterations_complete = 0;
    if (print_nodes_bool) { print_nodes(); }

    // Initialization conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << duration_string(start, stop) << ")\n";
}
void Mesh::_initialize(double _initial_contact_angle)
{
    // Constants
    constexpr auto PI = 3.141593265358979323846;
    const double ρ = 1.0 / sin(_initial_contact_angle);  // Spherical radius
    const double Δz = -ρ * cos(_initial_contact_angle);  // Vertical shift
    double volume_factor = 1.1;  // The factor to purposely overscale the volume

    // Coordinates (Polar, Azimuth | Cartesian)
    double θ, Δθ, φ, Δφ;
    double z, r, x, y;

    // Create middle point (φ = 0)
    current_nodes[res2][res2] = Node(0.0, 0.0, ρ + Δz);

    // Defining azimuth
    Δφ = _initial_contact_angle / res2;
    φ = Δφ;

    // Loop over radial segments of the spherical cap surface
    for (int j = 1; j <= res2; j++)  // Radial incrementer
    {
        // Defining cylindrical coordinates
        r = ρ * sin(φ);
        z = ρ * cos(φ) + Δz;

        // Defining the polar angle
        Δθ = (PI / 2) / (2 * j);
        θ = 0;

        // Loop over polar rotations about the z-axis
        for (int i = 0; i < 2 * j; i++)  // Polar incrementer
        {
            // Defining cartesian coordinates
            x = r * cos(θ);
            y = r * sin(θ);

            // Initialize node position vectors
            current_nodes[res2 - j + i][res2 - j] = Node( x,  y, z);  // Bottom
            current_nodes[res2 + j][res2 - j + i] = Node(-y,  x, z);  // Right
            current_nodes[res2 + j - i][res2 + j] = Node(-x, -y, z);  // Top
            current_nodes[res2 - j][res2 + j - i] = Node( y, -x, z);  // Left

            // Increment polar angle
            θ += Δθ;
        }

        // Increment azimuthal angle
        φ += Δφ;
    }

    // Scale node heights for a mesh volume ~ v_factor * expected volume
    volume_factor *= (current_droplet.volume_ratio / _calculate_volume());
    for (int i = 0; i < res; i++)
        for (int j = 0; j < res; j++)
            current_nodes[i][j].z *= volume_factor;

    // Update mesh characteristics
    _update_pressure();  // Update volume called within _update_pressure()
}

/*
 * Iterates the mesh _steps steps toward final geometry.
 * @brief   Iterate mesh's node arrays.
 */
void Mesh::iterate(int _steps)
{
    // Initialization start
    std::cout << "Iterating Mesh... ";
    time start = hrc::now();

    // Run function
    _iterate(_steps);

    // Print if enabled
    if (print_nodes_bool) { print_nodes(); }

    // Initialization conclusion
    time stop = hrc::now();
    std::cout << "Complete! (" << duration_string(start, stop) << ")\n";
}
void Mesh::_iterate(int _steps)
{
    // Variables
    Node mean, diff;
    double surface_contact_angle;

    // Iterate the mesh toward the expected volume
    for (int λ = 0; λ < _steps; λ++)
    {
        // Swap nodes in memory
        swap_nodes = previous_nodes;
        previous_nodes = current_nodes;
        current_nodes = swap_nodes;

        // Create current nodes from previous nodes
        for (int i = 0; i < res; i++)
            for (int j = 0; j < res; j++)
            {
                // Boundary points
                if (i == 0 || i == res1 || j == 0 || j == res1)
                {
                    // Find the node's contact angle
                    surface_contact_angle = _contact_angle(i, j);

                    // If node is on the printed region(s) and under the printed slip angle
                    if (current_surface.slips_on_printed(previous_nodes[i][j], surface_contact_angle))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i][j + 1]) / 3;
                            diff = previous_nodes[i + 1][j + 1] + _vector_gradient(i, j) * (previous_nodes[i + 1][j + 1].z * current_surface.printed_invtan);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i + 1][j - 1] + _vector_gradient(i, j) * (previous_nodes[i + 1][j - 1].z * current_surface.printed_invtan);
                        }

                        // Top-right boundary point
                        else if (i == res1 && j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i - 1][j] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i - 1][j - 1] + _vector_gradient(i, j) * (previous_nodes[i - 1][j - 1].z * current_surface.printed_invtan);
                        }

                        // Bottom-right boundary point
                        else if (i == res1 && j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i - 1][j] + previous_nodes[i][j + 1]) / 3;
                            diff = previous_nodes[i - 1][j + 1] + _vector_gradient(i, j) * (previous_nodes[i - 1][j + 1].z * current_surface.printed_invtan);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i - 1][j]) / 3;
                            diff = previous_nodes[i][j + 1] + _vector_gradient(i, j) * (previous_nodes[i][j + 1].z * current_surface.printed_invtan);
                        }

                        // Right side
                        else if (i == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i][j + 1] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i - 1][j] + _vector_gradient(i, j) * (previous_nodes[i - 1][j].z * current_surface.printed_invtan);
                        }

                        // Top side
                        else if (j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i - 1][j]) / 3;
                            diff = previous_nodes[i][j - 1] + _vector_gradient(i, j) * (previous_nodes[i][j - 1].z * current_surface.printed_invtan);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i][j + 1] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i + 1][j] + _vector_gradient(i, j) * (previous_nodes[i + 1][j].z * current_surface.printed_invtan);
                        }

                        // Create new node
                        current_nodes[i][j] = mean * α + diff * β;
                        current_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If node is on the unprinted region(s) and under the unprinted region slip angle
                    else if (current_surface.slips_on_unprinted(previous_nodes[i][j], surface_contact_angle))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i][j + 1]) / 3;
                            diff = previous_nodes[i + 1][j + 1] + _vector_gradient(i, j) * (previous_nodes[i + 1][j + 1].z * current_surface.unprinted_invtan);
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i + 1][j - 1] + _vector_gradient(i, j) * (previous_nodes[i + 1][j - 1].z * current_surface.unprinted_invtan);
                        }

                        // Top-right boundary point
                        else if (i == res1 && j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i - 1][j] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i - 1][j - 1] + _vector_gradient(i, j) * (previous_nodes[i - 1][j - 1].z * current_surface.unprinted_invtan);
                        }

                        // Bottom-right boundary point
                        else if (i == res1 && j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i - 1][j] + previous_nodes[i][j + 1]) / 3;
                            diff = previous_nodes[i - 1][j + 1] + _vector_gradient(i, j) * (previous_nodes[i - 1][j + 1].z * current_surface.unprinted_invtan);
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i - 1][j]) / 3;
                            diff = previous_nodes[i][j + 1] + _vector_gradient(i, j) * (previous_nodes[i][j + 1].z * current_surface.unprinted_invtan);
                        }

                        // Right side
                        else if (i == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i][j + 1] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i - 1][j] + _vector_gradient(i, j) * (previous_nodes[i - 1][j].z * current_surface.unprinted_invtan);
                        }

                        // Top side
                        else if (j == res1)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i + 1][j] + previous_nodes[i - 1][j]) / 3;
                            diff = previous_nodes[i][j - 1] + _vector_gradient(i, j) * (previous_nodes[i][j - 1].z * current_surface.unprinted_invtan);
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (previous_nodes[i][j] + previous_nodes[i][j + 1] + previous_nodes[i][j - 1]) / 3;
                            diff = previous_nodes[i + 1][j] + _vector_gradient(i, j) * (previous_nodes[i + 1][j].z * current_surface.unprinted_invtan);
                        }

                        // Create new node
                        current_nodes[i][j] = mean * α + diff * β;
                        current_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If the node's contact angle is steeper than the slip angles (pinned node)
                    else 
                    {
                        current_nodes[i][j] = previous_nodes[i][j];
                    }
                }

                // Internal points
                else
                {
                    current_nodes[i][j] = previous_nodes[i][j] + _force_net_normal(i, j) + _force_net_tangential(i, j);
                }
            }

        // Update mesh characteristics
        _update_pressure();  // Update volume called within _update_pressure()
    }
}



/***********************************************************
**                Characteristic Functions                **
***********************************************************/

/*
 * Calculates and updates the current volume of the mesh's surface.
 * @brief	Mesh volume functions.
 */
double Mesh::_update_volume()
{
    volume = _calculate_volume();
    return volume;
}
double Mesh::_calculate_volume()
{
    // Variables
    double V = 0.0;
    double h, dA;
    Node v1, v2, v3, v4;

    // Sum all the volume segments
    for (int i = 0; i < res1; i++)
        for (int j = 0; j < res1; j++)
        {
            // Assign variables
            v1 = current_nodes[i + 1][  j  ] - current_nodes[  i  ][  j  ];
            v2 = current_nodes[  i  ][j + 1] - current_nodes[  i  ][  j  ];
            v3 = current_nodes[i + 1][j + 1] - current_nodes[  i  ][j + 1];
            v4 = current_nodes[i + 1][j + 1] - current_nodes[i + 1][  j  ];

            // Calculate characteristics
            h = current_nodes[i][j].z + current_nodes[i + 1][j].z + current_nodes[i][j + 1].z + current_nodes[i + 1][j + 1].z;
            dA = (v1.x) * (v2.y) - (v2.x) * (v1.y) + (v3.x) * (v4.y) - (v4.x) * (v3.y);
            
            // Add new volume to total volume
            V += h * dA;
        }

    // Return volume
    return V * 0.125;  // Divide by 8 (1/4 for height, 1/2 for area, since it takes means of the values)
}

/*
 * Calculates and updates the current pressure of the mesh's surface.
 * Also updates volume when ran since it depends on volume.
 * @brief	Mesh pressure functions.
 */
double Mesh::_update_pressure()
{
    pressure = _calculate_pressure();
    return pressure;
}
double Mesh::_calculate_pressure()
{
    // Updates volume (just in case)
    _update_volume();

    // Iterate Gamma factor
    gamma += δ * current_droplet.contact_radius_cubed * (current_droplet.volume_ratio - volume);

    // Calculate pressure terms
    double base = current_droplet.volume_ratio / (volume * volume);
    double exponent = (1. + 0.1 * (current_droplet.volume_ratio / volume + volume / current_droplet.volume_ratio - 2.0)) / 3.0;

    // Return pressure
    return exp(gamma) * (σ / µ) * pow(base, exponent);
}



/***********************************************************
**                Approximation Functions                 **
***********************************************************/

/*
 * Returns the surface contact angle at point (i,j).
 * @note    Node(i, j) is only valid for boundary points.
 * @brief	The contact angle of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The contact angle IN RADIANS.
 */
double Mesh::_contact_angle(int i, int j)
{
    // Variables
    Node p1, p2, p3, p4, p5, p6;
    double approx_z;


    // Bottom-left edge point
    if (i == 0 && j == 0)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i + 1][  j  ]; p3 = previous_nodes[  i  ][j + 1];
        p4 = previous_nodes[i + 1][j + 1]; p5 = previous_nodes[i + 2][j + 1]; p6 = previous_nodes[i + 1][j + 2];
    }

    // Top-left edge point
    else if (i == 0 && j == res1)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i + 1][  j  ]; p3 = previous_nodes[  i  ][j - 1];
        p4 = previous_nodes[i + 1][j - 1]; p5 = previous_nodes[i + 2][j - 1]; p6 = previous_nodes[i + 1][j - 2];
    }

    // Top-right edge point
    else if (i == res1 && j == res1)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i - 1][  j  ]; p3 = previous_nodes[  i  ][j - 1];
        p4 = previous_nodes[i - 1][j - 1]; p5 = previous_nodes[i - 2][j - 1]; p6 = previous_nodes[i - 1][j - 2];
    }

    // Bottom-right edge point
    else if (i == res1 && j == 0)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i - 1][  j  ]; p3 = previous_nodes[  i  ][j + 1];
        p4 = previous_nodes[i - 1][j + 1]; p5 = previous_nodes[i - 2][j + 1]; p6 = previous_nodes[i - 1][j + 2];
    }

    // Bottom edge
    else if (j == 0)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i - 1][  j  ]; p3 = previous_nodes[i + 1][  j  ];
        p4 = previous_nodes[  i  ][j + 1]; p5 = previous_nodes[i - 1][j + 1]; p6 = previous_nodes[i + 1][j + 1];
    }

    // Right edge
    else if (i == res1)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[  i  ][j - 1]; p3 = previous_nodes[  i  ][j + 1];
        p4 = previous_nodes[i - 1][  j  ]; p5 = previous_nodes[i - 1][j - 1]; p6 = previous_nodes[i - 1][j + 1];
    }

    // Top edge
    else if (j == res1)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[i - 1][  j  ]; p3 = previous_nodes[i + 1][  j  ];
        p4 = previous_nodes[  i  ][j - 1]; p5 = previous_nodes[i - 1][j - 1]; p6 = previous_nodes[i + 1][j - 1];
    }

    // Left edge
    else if (i == 0)
    {
        p1 = previous_nodes[  i  ][  j  ]; p2 = previous_nodes[  i  ][j - 1]; p3 = previous_nodes[  i  ][j + 1];
        p4 = previous_nodes[i + 1][  j  ]; p5 = previous_nodes[i + 1][j - 1]; p6 = previous_nodes[i + 1][j + 1];
    }


    // Take a step in the gradient direction
    int sign = (p1.x * p1.x + p1.y * p1.y < p4.x * p4.x + p4.y * p4.y) ? 1 : -1;
    Node new_node = p1 - _vector_gradient(i, j) * grad_step * sign;


    // Endless variables
    double diagonal_factor, scale_factor;
    double matrix[6][7] = {
        { p1.x * p1.x, p1.x * p1.y, p1.y * p1.y, p1.x, p1.y, 1.0, p1.z },
        { p2.x * p2.x, p2.x * p2.y, p2.y * p2.y, p2.x, p2.y, 1.0, p2.z },
        { p3.x * p3.x, p3.x * p3.y, p3.y * p3.y, p3.x, p3.y, 1.0, p3.z },
        { p4.x * p4.x, p4.x * p4.y, p4.y * p4.y, p4.x, p4.y, 1.0, p4.z },
        { p5.x * p5.x, p5.x * p5.y, p5.y * p5.y, p5.x, p5.y, 1.0, p5.z },
        { p6.x * p6.x, p6.x * p6.y, p6.y * p6.y, p6.x, p6.y, 1.0, p6.z },
    };

    // 3rd Quadrant Triangle
    for (int diagonal = 0; diagonal < 6; diagonal++)
    {
        diagonal_factor = matrix[diagonal][diagonal];
        for (int row = diagonal; row < 6; row++)
        {
            scale_factor = diagonal_factor / matrix[row][0];
            for (int col = diagonal; col < 7; col++)
            {
                matrix[row][col] *= scale_factor;
                matrix[row][col] -= matrix[diagonal][col];
            }
        }
    }

    // 1st Quadrant Triangle
    for (int diagonal = 6 - 1; diagonal >= 0; diagonal--)
    {
        diagonal_factor = matrix[diagonal][diagonal];
        for (int row = diagonal - 1; row >= 0; row--)
        {
            scale_factor = diagonal_factor / matrix[row][0];
            for (int col = 7 - 1; col >= 0; col--)
            {
                matrix[row][col] *= scale_factor;
                matrix[row][col] -= matrix[diagonal][col];
            }
        }
    }

    // Return approximate z value at Node (x,y)
    approx_z = 
        matrix[0][6] * new_node.x * new_node.x +
        matrix[1][6] * new_node.x * new_node.y +
        matrix[2][6] * new_node.y * new_node.y +
        matrix[3][6] * new_node.x +
        matrix[4][6] * new_node.y +
        matrix[5][6];

    // Return the inverse tangent of the slope of the change in the node position
    return atan2(approx_z, -grad_step * sign);
}



/***********************************************************
**                    Vector Functions                    **
***********************************************************/

/*
 * Returns the normalized gradient vector at point (i,j).
 * @note    Node(i, j) is only valid for boundary points.
 * @brief	The normalized gradient vector of a node.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The normalized gradient.
 */
Node Mesh::_vector_gradient(int i, int j)
{
    // Variables
    const int n = res1;
    Node v;

    // Choose vector to rotate by 90 degrees
         if (i == 0 && j == 0) { v = previous_nodes[i + 1][j] - previous_nodes[i][j + 1]; }  // Bottom-left edge point
    else if (i == 0 && j == n) { v = previous_nodes[i][j - 1] - previous_nodes[i + 1][j]; }  // Top-left edge point
    else if (i == n && j == n) { v = previous_nodes[i - 1][j] - previous_nodes[i][j - 1]; }  // Top-right edge point
    else if (i == n && j == 0) { v = previous_nodes[i][j + 1] - previous_nodes[i - 1][j]; }  // Bottom-right edge point
    else if (j == 0)           { v = previous_nodes[i + 1][j] - previous_nodes[i - 1][j]; }  // Bottom edge
    else if (i == n)           { v = previous_nodes[i][j + 1] - previous_nodes[i][j - 1]; }  // Right edge
    else if (j == n)           { v = previous_nodes[i - 1][j] - previous_nodes[i + 1][j]; }  // Top edge
    else if (i == 0)           { v = previous_nodes[i][j - 1] - previous_nodes[i][j + 1]; }  // Left edge

    // Return normalized gradient vector
    return Node(v.y, -v.x, 0).normalize();
}

/*
 * Returns the normal vector at point(i, j) of the surface.
 * @note    Node(i, j) is only valid for non-boundary points.
 * @brief	Normal vector at surface.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The normal vector at that point.
 */
Node Mesh::_vector_normal(int i, int j)
{
    Node v1 = previous_nodes[i + 1][j] - previous_nodes[i - 1][j];
    Node v2 = previous_nodes[i][j + 1] - previous_nodes[i][j - 1];
    return cross_product(v1, v2).normalize();
}

/*
 * Returns the tangential vector to the surface centered at point (i,j).
 * @brief	Mesh tangential vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The tangent vector at that point.
 */
Node Mesh::_vector_tangent(int i, int j)
{
    Node v_cross = previous_nodes[i][j - 1] + previous_nodes[i][j + 1] + previous_nodes[i + 1][j] + previous_nodes[i - 1][j] - previous_nodes[i][j] * 4;
    return v_cross - v_cross.project(_vector_normal(i, j));
}

/*
 * Returns the mean curvature at point (i,j) by approximating second derivative.
 * @brief	Mean curvature vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The mean curvature vector at that point.
 */
Node Mesh::_vector_mean_curvature(int i, int j)
{
    // Variables
    Node v_i = previous_nodes[i + 1][j] - previous_nodes[i - 1][j];  // 2nd Difference in X
    Node v_j = previous_nodes[i][j + 1] - previous_nodes[i][j - 1];  // 2nd Difference in Y

    Node v_b = (previous_nodes[i][j - 1] - previous_nodes[i][j]).normalize();  // Bottom
    Node v_l = (previous_nodes[i - 1][j] - previous_nodes[i][j]).normalize();  // Left
    Node v_r = (previous_nodes[i + 1][j] - previous_nodes[i][j]).normalize();  // Right
    Node v_t = (previous_nodes[i][j + 1] - previous_nodes[i][j]).normalize();  // Left

    Node s1 = (v_l + v_r) * v_j.magnetude();
    Node s2 = (v_t + v_b) * v_i.magnetude();

    // Return mean curvature
    return (s1 + s2) * 0.5;
}



/***********************************************************
**                 Force Vector Functions                 **
***********************************************************/

/*
 * Returns the mesh curvature normal vector force at point (i,j).
 * @brief	Mesh curvature force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The curvature force at that point.
 */
Node Mesh::_force_curvature(int i, int j)
{
    Node normal_vector = _vector_normal(i, j);
    return _vector_mean_curvature(i, j).project(normal_vector) * (σ / µ);
}

/*
 * Returns the pressure normal force at point (i,j).
 * @brief	Mesh pressure force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The pressure force vector at that point.
 */
Node Mesh::_force_pressure(int i, int j)
{
    // Variables
    Node v1 = previous_nodes[i - 1][j] - previous_nodes[i][j - 1];
    Node v2 = previous_nodes[i - 1][j] - previous_nodes[i][j + 1];
    Node v3 = previous_nodes[i + 1][j] - previous_nodes[i][j - 1];
    Node v4 = previous_nodes[i + 1][j] - previous_nodes[i][j + 1];

    // Calculated cross products
    Node vec1 = cross_product(v1, v2);
    Node vec2 = cross_product(v3, v4);

    // Return normal pressure force
    return _vector_normal(i, j) * (pressure * 0.25) * (vec1.magnetude() + vec2.magnetude());
}

/*
 * Returns the tangential force to the surface centered at point (i,j).
 * @brief	Mesh tangential force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The tangent force at that point.
 */
Node Mesh::_force_tangential(int i, int j)
{
    return _vector_tangent(i, j) * (τ / µ);
}

/*
 * Returns the net normal force to the surface centered at point (i,j).
 * @brief	Mesh net normal force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The net normal force at that point.
 */
Node Mesh::_force_net_normal(int i, int j)
{
    return _force_pressure(i, j) + _force_curvature(i, j);
}

/*
 * Returns the net tangential force to the surface centered at point (i,j).
 * @brief	Mesh net tangential force vector.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	norm	Node	The net tangent force at that point.
 */
Node Mesh::_force_net_tangential(int i, int j)
{
    return _force_tangential(i, j);
}