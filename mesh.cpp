//
// File: mesh.cpp
// Author: Stanley Goodwin
// Creation Date: 6/16/2022
// Last Modified: 6/23/2022
// Credit to Kara Maki for skeleton code.
//
#include <iostream>
#include <fstream>
#include "mesh.h"
#include "misc.h"

#define PI 3.141593265358979323846
#define debug_printing true

using namespace std::chrono; // Ease of use



// Droplet variables (Remove after non-dimensionalization)
#define m_dV 3.0E-9
#define m_dR 1.9407025E-3
#define m_dκ (m_dV / (m_dR * m_dR * m_dR))



/******************************************************
**                Simulation Functions               **
**             Iteration & Initialization            **
******************************************************/

// Returns whether a point (i,j) is on the printed region
bool Mesh::OnPrintedRegion(int i, int j, int λ)
{
    double abs_y_val = m_dR * abs(_node_array[i][j][λ].y);

    bool r1 = (abs_y_val <= 0.25 * 0.001);
    bool r2 = (abs_y_val <= 1.25 * 0.001 && abs_y_val >= 0.75 * 0.001);

    return r1 || r2;
}

// Initializes the mesh nodes
void Mesh::InitializeNodes()
{
    // Function initialization
    std::cout << "Calculating Initial Node Mesh... ";
    auto start = high_resolution_clock::now();


    // Function constants
    const double θ = PI / 4;
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
        // Angle values [φ resets to 0 each radius increment]
        φ = 0.0;
        Δφ = (PI / 2) / (_res1 - 2 * j);

        for (int i = j; i <= _res1 - j; i++)  // Current angle index
        {
            // Node polar -> cartesian position calculation
            X = r * cos(φ);
            Y = r * sin(φ);
            Z = (j != 0) ? sqrt(k - (r * r)) - c : 0;

            // Initialize node position vectors
            _node_array[    i    ][    j    ][0] = Node( X,  Y, Z); // Bottom
            _node_array[_res1 - j][    i    ][0] = Node(-Y,  X, Z); // Right
            _node_array[_res1 - i][_res1 - j][0] = Node(-X, -Y, Z); // Top
            _node_array[    j    ][_res1 - i][0] = Node( Y, -X, Z); // Right

            // Increment angle
            φ += Δφ;
        }

        // Decrement radius
        r += Δr;
    }


    // Calculate mesh pressure
    _pressure = Pressure(0);

    // Function conclusion
    printTimeElapsed(start);
}

// Iterates the mesh
void Mesh::Iterate()
{
    // Function initialization
    std::cout << "Iterating Mesh... ";
    auto start = high_resolution_clock::now();


    // Function constants
    const double θ_c = 75 * PI / 180;  // Contact angle

    // Function Variables
    Node mean;  // The mean value of nodes
    Node diff;  // The change in the nodes


    // Iterate mesh toward final geometry

    for (int λ = 1; λ <= _total_iterations; λ++)  // Current iteration number
    {
        // Show current iteration percentage
        #if debug_printing
        printProgress((double)λ / _total_iterations);
        #endif

        // Create current nodes from previous nodes
        for (int i = 0; i < _resolution; i++)  // Current angle index
        {
            for (int j = 0; j < _resolution; j++)  // Current radius index
            {
                // Boundary points
                if (i == 0 || i == _res1 || j == 0 || j == _res1) {

                    // Bottom left point
                    if (i == 0 && j == 0)
                    {
                        _node_array[i][j][λ] = (_node_array[i][j + 1][λ - 1] + _node_array[i + 1][j][λ - 1]) / 2;
                        _node_array[i][j][λ].z = 0;  // Boundary condition
                    }

                    // Top left point
                    else if (i == 0 && j == _res1)
                    {
                        _node_array[i][j][λ] = (_node_array[i][j - 1][λ - 1] + _node_array[i + 1][j][λ - 1]) / 2;
                        _node_array[i][j][λ].z = 0;  // Boundary condition
                    }

                    // Top right point
                    else if (i == _res1 && j == _res1)
                    {
                        _node_array[i][j][λ] = (_node_array[i][j - 1][λ - 1] + _node_array[i - 1][j][λ - 1]) / 2;
                        _node_array[i][j][λ].z = 0;  // Boundary condition
                    }

                    // Bottom right point
                    else if (i == _res1 && j == 0)
                    {
                        _node_array[i][j][λ] = (_node_array[i][j + 1][λ - 1] + _node_array[i - 1][j][λ - 1]) / 2;
                        _node_array[i][j][λ].z = 0;  // Boundary condition
                    }

                    // If node is on the printed region(s)
                    else if (OnPrintedRegion(i, j, λ - 1))
                    {
                        // Bottom side
                        if (j == 0)
                        {
                            mean = (_node_array[i + 1][j][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i - 1][j][λ - 1]) / 3;
                            diff = _node_array[i][j + 1][λ - 1] - NormalVectorBottom(i, λ - 1) * (_node_array[i][j + 1][λ - 1].z / tan(θ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Right side
                        else if (i == _res1)
                        {
                            mean = (_node_array[i][j + 1][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i][j - 1][λ - 1]) / 3;
                            diff = _node_array[i - 1][j][λ - 1] - NormalVectorRight(j, λ - 1) * (_node_array[i - 1][j][λ - 1].z / tan(θ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Top side
                        else if (j == _res1)
                        {
                            mean = (_node_array[i + 1][j][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i - 1][j][λ - 1]) / 3;
                            diff = _node_array[i][j - 1][λ - 1] - NormalVectorTop(i, λ - 1) * (_node_array[i][j - 1][λ - 1].z / tan(θ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (_node_array[i][j + 1][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i][j - 1][λ - 1]) / 3;
                            diff = _node_array[i + 1][j][λ - 1] - NormalVectorLeft(j, λ - 1) * (_node_array[i + 1][j][λ - 1].z / tan(θ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }
                    }

                    // If node is not on the printed region(s) (pinned)
                    else {
                        _node_array[i][j][λ] = _node_array[i][j][λ - 1];
                    }
                }

                // Internal points
                else
                {
                    _node_array[i][j][λ] = _node_array[i][j][λ - 1] + NetNormalForce(i, j, λ - 1) + NetTangentialForce(i, j, λ - 1);
                }
            }
        }

        // Calculate mesh pressure
        _pressure = Pressure(λ);
    }


    // Function conclusion
    printTimeElapsed(start);
}



/*******************************************************
**                   Normal Vectors                   **
**       All normal vectors are of unit length.       **
*******************************************************/

Node Mesh::NormalVectorBottom(int i, int λ)
{
    int j = 0;
    Node diff2 = _node_array[i + 1][j][λ] - _node_array[i - 1][j][λ];
    return Node(-diff2.y, diff2.x, 0).normalize();
}
Node Mesh::NormalVectorRight(int j, int λ)
{
    int i = _res1;
    Node diff2 = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    return Node(-diff2.y, diff2.x, 0).normalize();
}
Node Mesh::NormalVectorTop(int i, int λ)
{
    int j = _res1;
    Node diff2 = _node_array[i + 1][j][λ] - _node_array[i - 1][j][λ];
    return Node(diff2.y, -diff2.x, 0).normalize();
}
Node Mesh::NormalVectorLeft(int j, int λ)
{
    int i = 0;
    Node diff2 = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    return Node(diff2.y, -diff2.x, 0).normalize();
}
Node Mesh::NormalVector(int i, int j, int λ)
{
    Node v1 = _node_array[i + 1][j][λ] - _node_array[i - 1][j][λ];
    Node v2 = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    return cross_product(v1, v2).normalize();
}



/***********************************************************
**                 Calculation Functions                  **
***********************************************************/

// Current Mesh Volume
double Mesh::Volume(int λ)
{
    // Initial volume
    double volume = 0;

    // Sum all the volume segments
    for (int i = 0; i < _res1; i++)
    {
        for (int j = 0; j < _res1; j++)
        {
            Node v1 = _node_array[i + 1][j][λ] - _node_array[i][j][λ];
            Node v2 = _node_array[i][j + 1][λ] - _node_array[i][j][λ];
            Node v3 = _node_array[i + 1][j + 1][λ] - _node_array[i][j + 1][λ];
            Node v4 = _node_array[i + 1][j + 1][λ] - _node_array[i + 1][j][λ];
            Node height = _node_array[i][j][λ] + _node_array[i + 1][j][λ] + _node_array[i][j + 1][λ] + _node_array[i + 1][j + 1][λ];

            volume += (
                (v1.x) * (v2.y) - (v2.x) * (v1.y) +
                (v3.x) * (v4.y) - (v4.x) * (v3.y)
            ) * height.z;
        }
    }

    // Return volume
    return volume / 8;
}

// Current Pressure
double Mesh::Pressure(int λ)
{
    // Current volume
    double volume = Volume(λ);

    // Iterate Gamma factor
    _Γ += δ * (m_dκ - volume) * m_dR * m_dR * m_dR;

    // Calculate pressure
    double base = m_dκ / (volume * volume);
    double expo = (1. + .1 * (m_dκ / volume + volume / m_dκ - 2)) / 3.;
    double pressure = exp(_Γ) * (σ / µ) * pow(base, expo);

    // Return pressure
    return pressure;
}


// Mesh Tangential Vector
Node Mesh::TangentPart(int i, int j, int λ)
{
    Node v_cross = _node_array[i][j - 1][λ] + _node_array[i][j + 1][λ] + _node_array[i + 1][j][λ] + _node_array[i - 1][j][λ] - _node_array[i][j][λ] * 4;
    Node v_normal = NormalVector(i, j, λ);
    return v_cross - v_cross.proj(v_normal);
}

// Mesh Tangential Force
Node Mesh::TangentialForce(int i, int j, int λ)
{
    return TangentPart(i, j, λ) * (τ / µ);
}

// Mesh Net Tangential Force
Node Mesh::NetTangentialForce(int i, int j, int λ)
{
    return TangentialForce(i, j, λ);
}


// Mesh Curvature Vector
Node Mesh::MeanCurvatureIntegral(int i, int j, int λ)
{
    Node v_i = _node_array[i + 1][j][λ] - _node_array[i - 1][j][λ];
    Node v_j = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    
    Node v_left  = (_node_array[i - 1][j][λ] - _node_array[i][j][λ]).normalize();
    Node v_right = (_node_array[i + 1][j][λ] - _node_array[i][j][λ]).normalize();
    Node v_up    = (_node_array[i][j + 1][λ] - _node_array[i][j][λ]).normalize();
    Node v_down  = (_node_array[i][j - 1][λ] - _node_array[i][j][λ]).normalize();

    Node s1 = (v_left + v_right) * v_j.det();
    Node s2 = (v_up + v_down) * v_i.det();
    
    return (s1 + s2) / 2;
}

// Mesh Curvature Force
Node Mesh::CurvatureForce(int i, int j, int λ)
{
    Node normal_vector = NormalVector(i, j, λ);
    return MeanCurvatureIntegral(i, j, λ).proj(normal_vector) * (σ / µ);
}

// Mesh Pressure Force
Node Mesh::PressureForce(int i, int j, int λ)
{
    Node v1 = _node_array[i - 1][j][λ] - _node_array[i][j - 1][λ];
    Node v2 = _node_array[i - 1][j][λ] - _node_array[i][j + 1][λ];
    Node v3 = _node_array[i + 1][j][λ] - _node_array[i][j - 1][λ];
    Node v4 = _node_array[i + 1][j][λ] - _node_array[i][j + 1][λ];

    Node vector1 = cross_product(v1, v2);
    Node vector2 = cross_product(v3, v4);

    double det1 = sqrt(dot_product(vector1, vector1));
    double det2 = sqrt(dot_product(vector2, vector2));

    double coeff = (_pressure / 4) * (det1 + det2);

    return NormalVector(i, j, λ) * coeff;
}

// Mesh Net Normal Force [Finished]
Node Mesh::NetNormalForce(int i, int j, int λ)
{
    return PressureForce(i, j, λ) + CurvatureForce(i, j, λ);
}



/******************************************************
**                 Misc Functions                    **
**        The remainder of functions of mesh         **
******************************************************/

// Print current droplet's nodes to file [TODO]
void Mesh::vprintCurrentMassInformation(std::string s, Droplet droplet)
{
    // drop variables
    std::string end = ".txt";
    std::ofstream myfile;
    std::string beg = "";
    beg = s;

    for (int count = 0; count <= _total_iterations; count += 1200) //for (int count = 0; count < (_total_iterations + 1); count = count + 100)
    {
        beg = s;
        myfile.open(beg.append("X").append(std::to_string(count)).append(end));

        for (int i = 0; i < _resolution; i++)
        {
            for (int j = 0; j < _resolution; j++)
            {
                if (j == _resolution - 1)
                    myfile << _node_array[i][j][count].x * m_dR << "\n";
                else
                    myfile << _node_array[i][j][count].x * m_dR << " ";
            }
        }

        myfile.close();
        beg = s;
        myfile.open(beg.append("Y").append(std::to_string(count)).append(end));


        for (int i = 0; i < _resolution; i++)
        {
            for (int j = 0; j < _resolution; j++)
            {
                if (j == _resolution - 1)
                    myfile << _node_array[i][j][count].y * m_dR << "\n";
                else
                    myfile << _node_array[i][j][count].y * m_dR << " ";
            }
        }

        myfile.close();
        beg = s;
        myfile.open(beg.append("Z").append(std::to_string(count)).append(end));

        for (int i = 0; i < _resolution; i++)
        {
            for (int j = 0; j < _resolution; j++)
            {
                if (j == _resolution - 1)
                    myfile << _node_array[i][j][count].z * m_dR << "\n";
                else
                    myfile << _node_array[i][j][count].z * m_dR << " ";
            }
        }

        myfile.close();
    }

}

// Mesh constructor
Mesh::Mesh()
{
    // Allocate memory for the node arrays
    for (int i = 0; i < _resolution; i++)
    {
        for (int j = 0; j < _resolution; j++)
        {
            _node_array[i][j] = new Node[_total_iterations + 1];
        }
    }
}
