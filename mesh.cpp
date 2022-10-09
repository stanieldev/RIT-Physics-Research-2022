//
// File: mesh.cpp
// Author: Stanley Goodwin
// Creation Date: 6/16/2022
// Last Modified: 7/4/2022
// Credit to Kara Maki for skeleton code.
//
#include <iostream>
#include <fstream>
#include <string>
#include "mesh.h"
#include "misc.h"

#define debug_printing true



/***********************************************************
**                 Simulation Functions                   **
***********************************************************/

// Returns whether a point (i,j) is on the printed region
bool Mesh::OnPrintedRegion(int i, int j)
{
    // Gap width parameters
    double w_g = 0.5 * 0.001 / _droplet.radius;
    double w_p = 0.5 * 0.001 / _droplet.radius;

    // Absolute positions
    double abs_y_val = abs(_prev_nodes[i][j].y);

    // Calculated parameters
    double half_p = w_p / 2;

    // Region booleans
    bool r1 = (abs_y_val <= half_p);
    bool r2 = (abs_y_val <= w_g + w_p + half_p && abs_y_val >= w_g + half_p);
    //bool r3 = (abs_y_val <= 2 * w_g + 2 * w_p + half_p && abs_y_val >= 2 * w_g + w_p + half_p);
    //bool r4 = (abs_y_val <= 3 * w_g + 3 * w_p + half_p && abs_y_val >= 3 * w_g + 2 * w_p + half_p);

    // Return if point is on printed region
    return r1 || r2;
}

// Initializes the mesh nodes
void Mesh::InitializeNodes()
{
    // Initialization
    std::cout << "Calculating Initial Node Mesh... ";
    auto start = start_timer;


    // Constants
    const double k = 1 / (sin(θ) * sin(θ));
    const double c = 1 / tan(θ);

    // Polar coordinates
    double r, Δr;  // Droplet radius & radial increment
    double φ, Δφ;  // Droplet angle & angular increment

    // Cartesian coordinates
    double X;  // Calculated value (cosine term)
    double Y;  // Calculated value (sine term)
    double Z;  // Calculated value (sqrt term)


    // Radius values
    r = 1.0;
    Δr = -2.0 / _res1;

    // Calculate initial mesh geometry
    for (int j = 0; j <= _res1 / 2; j++)  // Current radius index
    {
        // Angle values [φ resets to 0 for each radius increment]
        φ = 0.0;
        Δφ = (PI / 2) / (_res1 - 2 * j);

        for (int i = j; i <= _res1 - j; i++)  // Current angle index
        {
            // Node polar -> cartesian position calculation
            X = r * cos(φ);
            Y = r * sin(φ);
            Z = (j != 0) ? sqrt(k - (r * r)) - c : 0;

            // Initialize node position vectors
            _curr_nodes[    i    ][    j    ] = Node( X,  Y, Z);  // Bottom
            _curr_nodes[_res1 - j][    i    ] = Node(-Y,  X, Z);  // Right
            _curr_nodes[_res1 - i][_res1 - j] = Node(-X, -Y, Z);  // Top
            _curr_nodes[    j    ][_res1 - i] = Node( Y, -X, Z);  // Left

            // Increment angle
            φ += Δφ;
        }

        // Decrement radius
        r += Δr;
    }


    // Calculate mesh characteristics
    _volume = Volume();
    _pressure = Pressure();

    // Save current nodes to file
    PrintCurrent(0);

    // Conclusion
    stop_timer(start);
}

// Iterates the mesh
void Mesh::Iterate()
{
    // Initialization
    std::cout << "Iterating Mesh... ";
    auto start = start_timer;


    // Variables
    Node mean;  // The mean value of nodes
    Node diff;  // The change in the nodes


    // Iterate mesh toward final geometry
    for (int λ = 1; λ <= _nλ; λ++)  // Current iteration number
    {
        // Swap current & previous node memory addresses
        _swap_nodes = _prev_nodes;
        _prev_nodes = _curr_nodes;
        _curr_nodes = _swap_nodes;


        // Create current nodes from previous nodes
        for (int i = 0; i < _res; i++)  // Current angle index
        {
            for (int j = 0; j < _res; j++)  // Current radius index
            {
                // Boundary points
                if (i == 0 || i == _res1 || j == 0 || j == _res1) {

                    // Bottom left point
                    if (i == 0 && j == 0)
                    {
                        _curr_nodes[i][j] = (_prev_nodes[i][j + 1] + _prev_nodes[i + 1][j]) / 2;
                        _curr_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // Top left point
                    else if (i == 0 && j == _res1)
                    {
                        _curr_nodes[i][j] = (_prev_nodes[i][j - 1] + _prev_nodes[i + 1][j]) / 2;
                        _curr_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // Top right point
                    else if (i == _res1 && j == _res1)
                    {
                        _curr_nodes[i][j] = (_prev_nodes[i][j - 1] + _prev_nodes[i - 1][j]) / 2;
                        _curr_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // Bottom right point
                    else if (i == _res1 && j == 0)
                    {
                        _curr_nodes[i][j] = (_prev_nodes[i][j + 1] + _prev_nodes[i - 1][j]) / 2;
                        _curr_nodes[i][j].z = 0;  // Boundary condition
                    }

                    // If node is on the printed region(s)
                    else if (OnPrintedRegion(i, j))
                    {
                        // Bottom side
                        if (j == 0)
                        {
                            mean = (_prev_nodes[i + 1][j] + _prev_nodes[i][j] + _prev_nodes[i - 1][j]) / 3;
                            diff = _prev_nodes[i][j + 1] + NormalVector(i, j) * (_prev_nodes[i][j + 1].z / tan(θ_c));

                            _curr_nodes[i][j] = mean * α + diff * β;
                            _curr_nodes[i][j].z = 0;  // Boundary condition
                        }

                        // Right side
                        else if (i == _res1)
                        {
                            mean = (_prev_nodes[i][j + 1] + _prev_nodes[i][j] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i - 1][j] + NormalVector(i, j) * (_prev_nodes[i - 1][j].z / tan(θ_c));

                            _curr_nodes[i][j] = mean * α + diff * β;
                            _curr_nodes[i][j].z = 0;  // Boundary condition
                        }

                        // Top side
                        else if (j == _res1)
                        {
                            mean = (_prev_nodes[i + 1][j] + _prev_nodes[i][j] + _prev_nodes[i - 1][j]) / 3;
                            diff = _prev_nodes[i][j - 1] + NormalVector(i, j) * (_prev_nodes[i][j - 1].z / tan(θ_c));

                            _curr_nodes[i][j] = mean * α + diff * β;
                            _curr_nodes[i][j].z = 0;  // Boundary condition
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (_prev_nodes[i][j + 1] + _prev_nodes[i][j] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i + 1][j] + NormalVector(i, j) * (_prev_nodes[i + 1][j].z / tan(θ_c));

                            _curr_nodes[i][j] = mean * α + diff * β;
                            _curr_nodes[i][j].z = 0;  // Boundary condition
                        }
                    }

                    // If node is not on the printed region(s) (pinned)
                    else {
                        _curr_nodes[i][j] = _prev_nodes[i][j];
                    }
                }

                // Internal points
                else
                {
                    _curr_nodes[i][j] = _prev_nodes[i][j] + NetNormalForce(i, j) + NetTangentialForce(i, j);
                }
            }

            // Show current iteration percentage
            #if debug_printing
            print_progress((double)λ / _nλ);
            #endif
        }


        // Calculate mesh characteristics
        _volume = Volume();
        _pressure = Pressure();

        // Save current nodes to file
        if (λ == 1200) {
            PrintCurrent(λ);
        }
    }


    // Conclusion
    stop_timer(start);
}



/***********************************************************
**                 Calculation Functions                  **
***********************************************************/

// Normal vector at surface
Node Mesh::NormalVector(int i, int j)
{
    // Variables
    Node v1 = _prev_nodes[i + 1][j] - _prev_nodes[i - 1][j];
    Node v2 = _prev_nodes[i][j + 1] - _prev_nodes[i][j - 1];

    // Bottom edge
    if (j == 0) {
        return Node(v1.y, -v1.x, 0).normalize();
    }

    // Right edge
    else if (i == _res1) {
        return Node(v2.y, -v2.x, 0).normalize();
    }

    // Top edge
    else if (j == _res1) {
        return Node(-v1.y, v1.x, 0).normalize();
    }

    // Left edge
    else if (i == 0) {
        return Node(-v2.y, v2.x, 0).normalize();
    }

    // Internal nodes
    else {
        return cross_product(v1, v2).normalize();
    }
}

// Current Mesh Volume
double Mesh::Volume()
{
    // Variables
    double volume = 0;
    double area = 0;
    double height = 0;
    Node v1, v2, v3, v4;

    // Sum all the volume segments
    for (int i = 0; i < _res1; i++)
    {
        for (int j = 0; j < _res1; j++)
        {
            // Assign variables
            v1 = _curr_nodes[i + 1][  j  ] - _curr_nodes[  i  ][  j  ];
            v2 = _curr_nodes[  i  ][j + 1] - _curr_nodes[  i  ][  j  ];
            v3 = _curr_nodes[i + 1][j + 1] - _curr_nodes[  i  ][j + 1];
            v4 = _curr_nodes[i + 1][j + 1] - _curr_nodes[i + 1][  j  ];

            // Calculate lengths
            area = (v1.x) * (v2.y) - (v2.x) * (v1.y) + (v3.x) * (v4.y) - (v4.x) * (v3.y);
            height = (_curr_nodes[i][j] + _curr_nodes[i + 1][j] + _curr_nodes[i][j + 1] + _curr_nodes[i + 1][j + 1]).z;

            // Add new volume to total volume
            volume += area * height;
        }
    }

    // Return volume
    return volume / 8;
}

// Current Pressure
double Mesh::Pressure()
{
    // Variables
    double m_dκ = _droplet.κ;
    double volume = Volume();
    double pressure;
    double base;
    double expo;

    // Iterate Gamma factor
    _Γ += δ * (m_dκ - volume);

    // Calculate pressure
    base = m_dκ / (volume * volume);
    expo = (1. + .1 * (m_dκ / volume + volume / m_dκ - 2)) / 3.;
    pressure = exp(_Γ) * (σ / µ) * pow(base, expo);

    // Return pressure
    return pressure;
}



/***********************************************************
**                     Force Vectors                      **
***********************************************************/

// Mesh Tangential Vector
Node Mesh::TangentPart(int i, int j)
{
    // Variables
    Node v_cross = _prev_nodes[i][j - 1] + _prev_nodes[i][j + 1] + _prev_nodes[i + 1][j] + _prev_nodes[i - 1][j] - _prev_nodes[i][j] * 4;
    Node v_normal = NormalVector(i, j);

    // Return vector
    return v_cross - v_cross.proj(v_normal);
}

// Mesh Tangential Force
Node Mesh::TangentialForce(int i, int j)
{
    return TangentPart(i, j) * (τ / µ);
}

// Mesh Net Tangential Force
Node Mesh::NetTangentialForce(int i, int j)
{
    return TangentialForce(i, j);
}

// Mesh Curvature Vector
Node Mesh::MeanCurvatureIntegral(int i, int j)
{
    // Variables
    Node v_i = _prev_nodes[i + 1][j] - _prev_nodes[i - 1][j];  // 2nd Difference in X
    Node v_j = _prev_nodes[i][j + 1] - _prev_nodes[i][j - 1];  // 2nd Difference in Y

    Node v_b = (_prev_nodes[i][j - 1] - _prev_nodes[i][j]).normalize();  // Bottom
    Node v_l = (_prev_nodes[i - 1][j] - _prev_nodes[i][j]).normalize();  // Left
    Node v_r = (_prev_nodes[i + 1][j] - _prev_nodes[i][j]).normalize();  // Right
    Node v_t = (_prev_nodes[i][j + 1] - _prev_nodes[i][j]).normalize();  // Left

    Node s1 = (v_l + v_r) * v_j.det();
    Node s2 = (v_t + v_b) * v_i.det();
    
    // Return mean curvature
    return (s1 + s2) / 2;  // Second derivative approximation
}

// Mesh Curvature Force
Node Mesh::CurvatureForce(int i, int j)
{
    // Variables
    Node normal_vector = NormalVector(i, j);

    // Return mean curvature normal vector
    return MeanCurvatureIntegral(i, j).proj(normal_vector) * (σ / µ);
}

// Mesh Pressure Force
Node Mesh::PressureForce(int i, int j)
{
    // Variables
    Node v1 = _prev_nodes[i - 1][j] - _prev_nodes[i][j - 1];
    Node v2 = _prev_nodes[i - 1][j] - _prev_nodes[i][j + 1];
    Node v3 = _prev_nodes[i + 1][j] - _prev_nodes[i][j - 1];
    Node v4 = _prev_nodes[i + 1][j] - _prev_nodes[i][j + 1];

    // Calculated cross products
    Node vec1 = cross_product(v1, v2);
    Node vec2 = cross_product(v3, v4);

    // Return normal pressure force
    return NormalVector(i, j) * (_pressure / 4) * (vec1.det() + vec2.det());
}

// Mesh Net Normal Force [Finished]
Node Mesh::NetNormalForce(int i, int j)
{
    return PressureForce(i, j) + CurvatureForce(i, j);
}



/***********************************************************
**                     Misc Functions                     **
***********************************************************/

// Printing current nodes to output file
void Mesh::PrintCurrent(int λ)
{
    // File variables
    std::ofstream fileX, fileY, fileZ;
    std::string directory = "C:\\Users\\sfg99\\Code\\Summer Research\\Matlab\\data_generated\\";
    std::string file_name = "Xitwocompleted41";
    std::string file_end = ".txt";
    
    // Open file streams
    fileX.open(directory + file_name + "X" + std::to_string(λ) + file_end);
    fileY.open(directory + file_name + "Y" + std::to_string(λ) + file_end);
    fileZ.open(directory + file_name + "Z" + std::to_string(λ) + file_end);
    
    // Print all nodes to separate files
    Node current_node;
    for (int i = 0; i < _res; i++)
    {
        for (int j = 0; j < _res; j++)
        {
            // Set current node
            current_node = _curr_nodes[i][j] * _droplet.radius;

            // Send data to the files
            fileX << current_node.x << " ";
            fileY << current_node.y << " ";
            fileZ << current_node.z << " ";
        }

        // Add new line
        fileX << "\n";
        fileY << "\n";
        fileZ << "\n";
    }

    // Close file streams
    fileX.close();
    fileY.close();
    fileZ.close();
}

// Mesh constructor
Mesh::Mesh(Droplet droplet, int resolution, int total_iteration_count)
{
    // User parameters
    _droplet = droplet;
    _nλ = total_iteration_count;
    _res = resolution;
    _res1 = _res - 1;

    // Heap allocation
    /*Node** _prev_array;
    _prev_array = new Node*[_res];
    for (int i = 0; i < _res; i++)
    {
        _prev_array[i] = new Node[_res];
    }*/

    // Extra calculations
    δ *= droplet.radius * droplet.radius * droplet.radius;
}
