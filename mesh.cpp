//
// File: mesh.cpp
// Author: Stanley Goodwin
// Creation Date: 6/16/2022
// Last Modified: 6/16/2022
// Credit to Kara Maki for skeleton code.
//
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <chrono>
#include "mesh.h"

#define PI 3.141593265358979323846
#define debug_printing true





// Droplet variables (Remove after non-dimensionalization)
#define m_dV 3.0E-9
#define m_dR 1.9407025E-3



// Progress bar function
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
void printProgress(double percentage) { // https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
    //percentage -= 1;
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}



void printTimeElapsed(std::chrono::high_resolution_clock::time_point start_time)
{
    using namespace std::chrono;

    // The time the function is called (effectively the end-time)
    auto stop_time = high_resolution_clock::now();

    // Casting and reducing until variables are made
    auto µ = duration_cast<microseconds>(stop_time - start_time);
    auto ms = duration_cast<milliseconds>(µ);
    µ -= duration_cast<microseconds>(ms);
    auto s = duration_cast<seconds>(ms);
    ms -= duration_cast<milliseconds>(s);
    auto m = duration_cast<minutes>(s);
    s -= duration_cast<seconds>(m);
    auto h = duration_cast<hours>(m);
    m -= duration_cast<minutes>(h);

    // Console output
    std::cout << "Complete! (";
    if (h.count() != 0) {
        std::cout << h.count() << "h, " << m.count() << "m";
    }
    else if (m.count() != 0) {
        std::cout << m.count() << "m, " << s.count() << "s";
    }
    else if (s.count() != 0) {
        std::cout << s.count() + ms.count() / 1000.0 << "s";
    }
    else if (ms.count() != 0) {
        std::cout << ms.count() + µ.count() / 1000.0 << "ms";
    }
    else {
        std::cout << µ.count() << "us";
    }
    std::cout << ")\n";
}



/******************************************************
**                Simulation Functions               **
**             Iteration & Initialization            **
******************************************************/



void Mesh::CInitNodes()
{
    // Function initialization
    std::cout << "Calculating Initial Node Mesh... ";
    auto start = std::chrono::high_resolution_clock::now();


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

            // TODO, REMOVE RADIUS FROM THESE EQUATIONS
            X *= m_dR; Y *= m_dR; Z *= m_dR;

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


    // Function conclusion
    printTimeElapsed(start);
}


void CBoundary(int λ);
void CInterior(int λ);









void Mesh::Iterate()
{
    // Function initialization
    std::cout << "Iterating Mesh... ";
    auto start = std::chrono::high_resolution_clock::now();


    // Function definitions
    const double φ_c = -75 * PI / 180;  // Contact angle

    int λ;  // Current iteration (pseudo-time)
    int i;  // Current angle index
    int j;  // Current radius index

    Node mean;  // The mean value of nodes
    Node diff;  // The change in the nodes


    // Iterate mesh toward final geometry
    for (λ = 1; λ <= _total_iterations; λ++)  // Mesh-Iteration iterator
    {
        // Show current iteration percentage
#if debug_printing
        printProgress((double)λ / _total_iterations);
#endif

        // Create current points from previous points
        for (i = 0; i < _resolution; i++)  // Angle iterator
        {
            for (j = 0; j < _resolution; j++)  // Radius iterator
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
                    else if (
                        abs(_node_array[i][j][λ - 1].y) <= 0.25 * 0.001 ||
                        (abs(_node_array[i][j][λ - 1].y) <= 1.25 * 0.001) && (abs(_node_array[i][j][λ - 1].y) >= 0.75 * 0.001)
                    ) {
                        // Bottom side
                        if (j == 0) {
                            mean = (_node_array[i + 1][j][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i - 1][j][λ - 1]) / 3;
                            diff = _node_array[i][j + 1][λ - 1] + NormalVectorBottom(i, λ - 1) * (_node_array[i][j + 1][λ - 1].z / tan(φ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Right side
                        else if (i == _res1)
                        {
                            mean = (_node_array[i][j + 1][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i][j - 1][λ - 1]) / 3;
                            diff = _node_array[i - 1][j][λ - 1] + NormalVectorRight(j, λ - 1) * (_node_array[i - 1][j][λ - 1].z / tan(φ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Top side
                        else if (j == _res1)
                        {
                            mean = (_node_array[i + 1][j][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i - 1][j][λ - 1]) / 3;
                            diff = _node_array[i][j - 1][λ - 1] + NormalVectorTop(i, λ - 1) * (_node_array[i][j - 1][λ - 1].z / tan(φ_c));

                            _node_array[i][j][λ] = mean * α + diff * β;
                            _node_array[i][j][λ].z = 0;  // Boundary condition
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (_node_array[i][j + 1][λ - 1] + _node_array[i][j][λ - 1] + _node_array[i][j - 1][λ - 1]) / 3;
                            diff = _node_array[i + 1][j][λ - 1] + NormalVectorLeft(j, λ - 1) * (_node_array[i + 1][j][λ - 1].z / tan(φ_c));

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
                    _node_array[i][j][λ].x = _node_array[i][j][λ - 1].x + (NetNormalForce(i, j, λ - 1).x + NetTangentialForce(i, j, λ - 1).x) * (Δt / µ);
                    _node_array[i][j][λ].y = _node_array[i][j][λ - 1].y + (NetNormalForce(i, j, λ - 1).y + NetTangentialForce(i, j, λ - 1).y) * (Δt / µ);
                    _node_array[i][j][λ].z = _node_array[i][j][λ - 1].z + (NetNormalForce(i, j, λ - 1).z + NetTangentialForce(i, j, λ - 1).z) * (Δt / µ);
                    //_node_array[i][j][λ] = _node_array[i][j][λ - 1]   + (NetNormalForce(i, j, λ - 1)   + NetTangentialForce(i, j, λ - 1)  ) * (Δt / µ);
                }
            }
        }

        vPressure(λ - 1);
    }


    // Function conclusion
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Complete! (" << duration.count() << "s)\n";
}







/******************************************************
**                 Vector Arithmetic                 **
**         The dot product and cross product         **
******************************************************/

double dot_product(Node v1, Node v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
Node cross_product(Node v1, Node v2)
{
    double i = v1.y * v2.z - v1.z * v2.y;
    double j = v1.z * v2.x - v1.x * v2.z;
    double k = v1.x * v2.y - v1.y * v2.x;
    return Node(i, j, k);
};



/*******************************************************
**                   Normal Vectors                   **
**  The normal vectors for all places on the surface  **
*******************************************************/






Node Mesh::NormalVectorBottom(int i, int λ)
{
    Node diff2 = _node_array[i + 1][0][λ] - _node_array[i - 1][0][λ];
    return Node(-diff2.y, diff2.x, 0).normalize();  // +90 Degree rotation
}
Node Mesh::NormalVectorRight(int j, int λ)
{
    Node diff2 = _node_array[_res1][j + 1][λ] - _node_array[_res1][j - 1][λ];
    return Node(-diff2.y, diff2.x, 0).normalize();  // +90 Degree rotation
}
Node Mesh::NormalVectorTop(int i, int λ)
{
    Node diff2 = _node_array[i + 1][_res1][λ] - _node_array[i - 1][_res1][λ];
    return Node(diff2.y, -diff2.x, 0).normalize();  // -90 Degree rotation
}
Node Mesh::NormalVectorLeft(int j, int λ)
{
    Node diff2 = _node_array[0][j + 1][λ] - _node_array[0][j - 1][λ];
    return Node(diff2.y, -diff2.x, 0).normalize();  // -90 Degree rotation
}
Node Mesh::NormalVector(int i, int j, int λ)
{
    Node v1 = _node_array[i + 1][j][λ] - _node_array[i - 1][j][λ];
    Node v2 = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    return cross_product(v1, v2).normalize();
}


/******************************************************
**              The Internal Conditions              **
**     The normal vectors for all each dimension     **
******************************************************/

double Mesh::vCurrentVolume(int λ)
{
    double v_current = 0;

    for (int i = 0; i < _res1; i++)
    {
        for (int j = 0; j < _res1; j++)
        {
            v_current += (
                (_node_array[i + 1][  j  ][λ].x - _node_array[  i  ][  j  ][λ].x) *
                (_node_array[  i  ][j + 1][λ].y - _node_array[  i  ][  j  ][λ].y) -
                (_node_array[  i  ][j + 1][λ].x - _node_array[  i  ][  j  ][λ].x) *
                (_node_array[i + 1][  j  ][λ].y - _node_array[  i  ][  j  ][λ].y) +
                (_node_array[i + 1][j + 1][λ].x - _node_array[  i  ][j + 1][λ].x) *
                (_node_array[i + 1][j + 1][λ].y - _node_array[i + 1][  j  ][λ].y) -
                (_node_array[i + 1][j + 1][λ].x - _node_array[i + 1][  j  ][λ].x) *
                (_node_array[i + 1][j + 1][λ].y - _node_array[  i  ][j + 1][λ].y)
            ) * (_node_array[  i  ][  j  ][λ].z + _node_array[i + 1][  j  ][λ].z + 
                 _node_array[  i  ][j + 1][λ].z + _node_array[i + 1][j + 1][λ].z
            );
        }
    }
    return v_current / 8;
}
double Mesh::vPressure(int λ)
{
    double pr;
    double delta_p;

    pr = m_PressureFactor;

    delta_p = δ * (m_dV - vCurrentVolume(λ));
    pr += delta_p;

    m_dP = exp(pr) * σ * pow(m_dV / (vCurrentVolume(λ) * vCurrentVolume(λ)),
        (1. + .1 * ((m_dV / vCurrentVolume(λ) + vCurrentVolume(λ) / m_dV) - 2)) / 3.);

    m_PressureFactor = pr;
    return m_dP;
}

Node Mesh::MeanCurvature(int i, int j, int λ)
{
    Node v1 = _node_array[i][j + 1][λ] - _node_array[i][j - 1][λ];
    Node v4 = _node_array[i - 1][j][λ] - _node_array[i + 1][j][λ];

    Node v2 = _node_array[i - 1][j][λ] - _node_array[i][j][λ];
    Node v3 = _node_array[i + 1][j][λ] - _node_array[i][j][λ];
    Node v5 = _node_array[i][j + 1][λ] - _node_array[i][j][λ];
    Node v6 = _node_array[i][j - 1][λ] - _node_array[i][j][λ];

    double m_left = sqrt(dot_product(v1, v1) / dot_product(v2, v2));
    double m_right = sqrt(dot_product(v1, v1) / dot_product(v3, v3));
    double m_up = sqrt(dot_product(v4, v4) / dot_product(v5, v5));
    double m_down = sqrt(dot_product(v4, v4) / dot_product(v6, v6));

    return (v2 * m_left + v3 * m_right + v5 * m_up + v6 * m_down) / 2;
}
double Mesh::PressureForce(int i, int j, int λ)
{
    /******************************************************
    **                      FORMULA                      **
    **   P/4 * {|(Vi,j-1 - Vi-1,j)X(Vi-1,j - Vi,j+1 )|   **
    **   + |(Vi,j-1 - Vi+1,j)X(Vi+1,j - Vi,j+1 )|}       **
    ******************************************************/

    double p;

    double vYZ = (_node_array[i][j - 1][λ].y - _node_array[i - 1][j][λ].y) * (_node_array[i - 1][j][λ].z - _node_array[i][j + 1][λ].z);
    double vZY = (_node_array[i][j - 1][λ].z - _node_array[i - 1][j][λ].z) * (_node_array[i - 1][j][λ].y - _node_array[i][j + 1][λ].y);
    double vXZ = (_node_array[i][j - 1][λ].x - _node_array[i - 1][j][λ].x) * (_node_array[i - 1][j][λ].z - _node_array[i][j + 1][λ].z);
    double vZX = (_node_array[i][j - 1][λ].z - _node_array[i - 1][j][λ].z) * (_node_array[i - 1][j][λ].x - _node_array[i][j + 1][λ].x);
    double vXY = (_node_array[i][j - 1][λ].x - _node_array[i - 1][j][λ].x) * (_node_array[i - 1][j][λ].y - _node_array[i][j + 1][λ].y);
    double vYX = (_node_array[i][j - 1][λ].y - _node_array[i - 1][j][λ].y) * (_node_array[i - 1][j][λ].x - _node_array[i][j + 1][λ].x);

    double vYZ2 = (_node_array[i][j - 1][λ].y - _node_array[i + 1][j][λ].y) * (_node_array[i + 1][j][λ].z - _node_array[i][j + 1][λ].z);
    double vZY2 = (_node_array[i][j - 1][λ].z - _node_array[i + 1][j][λ].z) * (_node_array[i + 1][j][λ].y - _node_array[i][j + 1][λ].y);
    double vXZ2 = (_node_array[i][j - 1][λ].x - _node_array[i + 1][j][λ].x) * (_node_array[i + 1][j][λ].z - _node_array[i][j + 1][λ].z);
    double vZX2 = (_node_array[i][j - 1][λ].z - _node_array[i + 1][j][λ].z) * (_node_array[i + 1][j][λ].x - _node_array[i][j + 1][λ].x);
    //double vXY2 = (_node_array[i][j - 1][λ].x - _node_array[i + 1][j][λ].x) * (_node_array[i + 1][j][λ].y - _node_array[i][j + 1][λ].y);
    //double vYX2 = (_node_array[i][j - 1][λ].y - _node_array[i + 1][j][λ].y) * (_node_array[i + 1][j][λ].x - _node_array[i][j + 1][λ].x);

    //Node v1 = _node_array[i][j - 1][λ] - _node_array[i - 1][j][λ];
    //Node v2 = _node_array[i + 1][j][λ] - _node_array[i][j + 1][λ];
    

    Node vector = Node(
        (vXY - vYX),
        (vYZ - vZY),
        (vXZ - vZX)
    );
    Node vector2 = Node(
        (vYZ2 - vZY2),
        (vXZ2 - vZX2),
        (vYZ2 - vZY2)
        //(vXY2 - vYX2)
    );

    p = vPressure(λ) / 4 * (
        sqrt(dot_product(vector, vector)) +
        sqrt(dot_product(vector2, vector2))
    );

    return p;
}
double Mesh::SigmaForce(int i, int j, int λ)
{
    return dot_product(NormalVector(i, j, λ), MeanCurvature(i, j, λ));
}
Node Mesh::NetNormalForce(int i, int j, int λ)
{
    double pressure_force = PressureForce(i, j, λ);
    double curvature_force = SigmaForce(i, j, λ);
    return NormalVector(i, j, λ) * (pressure_force + σ * curvature_force);
}

Node Mesh::TangentPart(int i, int j, int λ)
{
    Node v1 = _node_array[i][j - 1][λ] + _node_array[i][j + 1][λ] + _node_array[i + 1][j][λ] + _node_array[i - 1][j][λ] - _node_array[i][j][λ] * 4;
    Node v2 = NormalVector(i, j, λ);

    double dot_p = dot_product(v1, v2);

    return v1 - v2 * dot_p;
}
Node Mesh::NetTangentialForce(int i, int j, int λ)
{
    return TangentPart(i, j, λ) * τ;
}



/******************************************************
**                 Misc Functions                    **
**        The remainder of functions of Mesh         **
******************************************************/

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
                    myfile << _node_array[i][j][count].x << "\n";
                else
                    myfile << _node_array[i][j][count].x << " ";
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
                    myfile << _node_array[i][j][count].y << "\n";
                else
                    myfile << _node_array[i][j][count].y << " ";
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
                    myfile << _node_array[i][j][count].z << "\n";
                else
                    myfile << _node_array[i][j][count].z << " ";
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
            _normal_vector[i][j] = new Node[_total_iterations + 1];
        }
    }















    /* OLD STUFF TO REPLACE / REUSE */
    m_PressureFactor = 0;
}