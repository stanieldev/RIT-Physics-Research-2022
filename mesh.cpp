/**
 * @file	mesh.cpp
 * @brief	The mesh function definitions.
 *
 * @author	Stanley Goodwin
 *          Kara Maki       (previous versions)
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/16/2022
 * Last Modified: 7/14/2022
 */
#include <iostream>
#include <fstream>
#include <string>
#include "fmath.h"
#include "misc.h"
#include "mesh.h"

#define debug_printing false



/***********************************************************
**                  Simulation Functions                  **
***********************************************************/

/*
* Function description:
*   Initializes the mesh node array to prepare it for iteration.
*   Creates an initial mesh in the shape of a spherical cap.
*
* Resources:
*   https://en.wikipedia.org/wiki/Spherical_coordinate_system
* 
* Notes:
*   Currently only works with an odd number of points.
*/
void Mesh::initialize()
{
    // Initialization
    printf("Calculating Initial Node Mesh... ");
    time start = hrc::now();


    // Constants
    const int n = _res / 2;  // Resolution midpoint (rounds down)

    // Spherical coordinates
    double θ, Δθ;  // Polar angle
    double φ, Δφ;  // Azimuthal angle

    // Derived coordinates (spherical -> cartesian where ρ = 1)
    double z;  // The azimuth cosine component of the droplet radius
    double r;  // The azimuth sine component of the droplet radius
    double x;  // The polar cosine component of r
    double y;  // The polar sine component of r
    

    // Defining the azimuthal angle
    Δφ = _θi / n;
    φ = (_res % 2) ? Δφ : Δφ / 2;  // Angular offset

    // Create middle point (if resolution is odd)
    if (_res % 2) { _curr_nodes[n][n] = Node(0, 0, 1 / sin(_θi) - 1 / tan(_θi)); }
    
    // Loop over radial segments of the spherical cap surface
    for (int j = 0; j < n; j++)  // Radial incrementer
    {
        // Defining variables
        r = sin(φ) / sin(_θi);
        z = cos(φ) / sin(_θi) - 1 / tan(_θi);

        // Defining the polar angle
        Δθ = (PI / 2) / (2 * j + 2);
        θ = 0;

        // Loop over polar rotations about the z-axis
        for (int i = 0; i < 2 * j + 2; i++)  // Angular incrementer
        {
            // Defining variables
            x = r * cos(θ);
            y = r * sin(θ);

            // Initialize node position vectors
            _curr_nodes[n - j - 1 + i][n - j - 1    ] = Node( x,  y, z);  // Bottom
            _curr_nodes[n + j + 1    ][n - j - 1 + i] = Node(-y,  x, z);  // Right
            _curr_nodes[n + j + 1 - i][n + j + 1    ] = Node(-x, -y, z);  // Top
            _curr_nodes[n - j - 1    ][n + j + 1 - i] = Node( y, -x, z);  // Left

            // Increment polar angle
            θ += Δθ;
        }

        // Increment azimuthal angle
        φ += Δφ;
    }


    // Calculate mesh characteristics
    _volume = Volume();
    _pressure = Pressure();

    // Save current nodes to file
    print_current(0);

    // Conclusion
    print_duration(start);
}

/*
* Function description:
*   Iterates the mesh toward the final geometry of the droplet.
*
* Notes:
*   Currently only works with an odd number of points.
*   New contact angle stuff (TODO)
*/
void Mesh::iterate()
{
    // Initialization
    printf("Iterating Mesh... ");
    time start = hrc::now();

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

                    // If node is on the printed region(s)
                    //printf("%f, ", contact_angle(i, j));
                    if (contact_angle(i, j) < θ_c && OnPrintedRegion(i, j))
                    {
                        // Bottom-left boundary point
                        if (i == 0 && j == 0)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i + 1][j] + _prev_nodes[i][j + 1]) / 3;
                            diff = _prev_nodes[i + 1][j + 1] + vector_gradient(i, j) * (_prev_nodes[i + 1][j + 1].z / tan(θ_c));  // change top ThetaD for Non-printed
                        }

                        // Top-left boundary point
                        else if (i == 0 && j == _res1)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i + 1][j] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i + 1][j - 1] + vector_gradient(i, j) * (_prev_nodes[i + 1][j - 1].z / tan(θ_c));
                        }

                        // Top-right boundary point
                        else if (i == _res1 && j == _res1)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i - 1][j] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i - 1][j - 1] + vector_gradient(i, j) * (_prev_nodes[i - 1][j - 1].z / tan(θ_c));
                        }

                        // Bottom-right boundary point
                        else if (i == _res1 && j == 0)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i - 1][j] + _prev_nodes[i][j + 1]) / 3;
                            diff = _prev_nodes[i - 1][j + 1] + vector_gradient(i, j) * (_prev_nodes[i - 1][j + 1].z / tan(θ_c));
                        }

                        // Bottom side
                        else if (j == 0)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i + 1][j] + _prev_nodes[i - 1][j]) / 3;
                            diff = _prev_nodes[i][j + 1] + vector_gradient(i, j) * (_prev_nodes[i][j + 1].z / tan(θ_c));
                        }

                        // Right side
                        else if (i == _res1)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i][j + 1] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i - 1][j] + vector_gradient(i, j) * (_prev_nodes[i - 1][j].z / tan(θ_c));
                        }

                        // Top side
                        else if (j == _res1)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i + 1][j] + _prev_nodes[i - 1][j]) / 3;
                            diff = _prev_nodes[i][j - 1] + vector_gradient(i, j) * (_prev_nodes[i][j - 1].z / tan(θ_c));
                        }

                        // Left side
                        else if (i == 0)
                        {
                            mean = (_prev_nodes[i][j] + _prev_nodes[i][j + 1] + _prev_nodes[i][j - 1]) / 3;
                            diff = _prev_nodes[i + 1][j] + vector_gradient(i, j) * (_prev_nodes[i + 1][j].z / tan(θ_c));
                        }

                        // Create new node
                        _curr_nodes[i][j] = mean * α + diff * β;
                        _curr_nodes[i][j].z = 0;  // Boundary condition
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
        if (λ == _nλ) {
            print_current(λ);
        }
    }

    // Conclusion
    print_duration(start);
}

/*
* Function description:
*   Prints the contents of the node arrays to a text file.
*
* @param	λ	int	  Current iteration number.
* 
* Notes:
*   Currently only works with an odd number of points.
*   New contact angle stuff (TODO)
*/
void Mesh::print_current(int λ)
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
            current_node = _curr_nodes[i][j];  //  * _droplet.contact_radius

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



/***********************************************************
**                     Mesh Functions                     **
***********************************************************/

/**
 * Initializes a node using user parameters.
 *
 * @brief	Mesh constructor.
 * @param	droplet	            Droplet	 The droplet to simulate.
 * @param	resolution	            int	 Node resolution number.
 * @param	total_iteration_count	int	 Max iteration number.
 */
Mesh::Mesh()
{
    // Heap allocation
    /*Node** _prev_array;
    _prev_array = new Node*[_res];
    for (int i = 0; i < _res; i++)
    {
        _prev_array[i] = new Node[_res];
    }*/

    // Extra calculations
    δ *= _droplet.contact_radius * _droplet.contact_radius * _droplet.contact_radius;
}

/**
 * Returns the approximate volume of the surface at the time of execution.
 *
 * @brief	Current mesh volume.
 * @return	volume	double	The volume of the current surface.
 */
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

/**
 * Returns the approximate pressure of the surface at the time of execution.
 *
 * @brief	Current pressure.
 * @return	presure	double	The pressure of the current surface.
 */
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
 **                   Boundary Functions                   **
 ***********************************************************/

/**
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
    int n = _res / 2;
    Node v;

    // Bottom-left edge point
    if (i == 0 && j == 0)
    {
        v = _prev_nodes[i + 1][j] - _prev_nodes[i][j + 1];
    }

    // Top-left edge point
    else if (i == 0 && j == _res1)
    {
        v = _prev_nodes[i][j - 1] - _prev_nodes[i + 1][j];
    }

    // Top-right edge point
    else if (i == _res1 && j == _res1)
    {
        v = _prev_nodes[i - 1][j] - _prev_nodes[i][j - 1];
    }

    // Bottom-right edge point
    else if (i == _res1 && j == 0)
    {
        v = _prev_nodes[i][j + 1] - _prev_nodes[i - 1][j];
    }

    // Bottom edge
    if (j == 0) {
        v = _prev_nodes[i + 1][j] - _prev_nodes[i - 1][j];
    }

    // Right edge
    else if (i == _res1) {
        v = _prev_nodes[i][j + 1] - _prev_nodes[i][j - 1];
    }

    // Top edge
    else if (j == _res1) {
        v = _prev_nodes[i - 1][j] - _prev_nodes[i + 1][j];
    }

    // Left edge
    else if (i == 0) {
        v = _prev_nodes[i][j - 1] - _prev_nodes[i][j + 1];
    }

    // Return normalized gradient vector
    return Node(v.y, -v.x, 0).normalize();
}


double Mesh::contact_angle(int i, int j)
{
    // Variables
    Node p1, p2, p3, p4, p5, p6;

    // Bottom-left edge point
    if (i == 0 && j == 0)
    {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i + 1][j];
        p3 = _prev_nodes[i + 1][j + 1];
        p4 = _prev_nodes[i + 2][j];
        p5 = _prev_nodes[i + 2][j + 1];
        p6 = _prev_nodes[i + 2][j + 2];

        printf(
            "(%f, %f, %f), \n",
            p1.x, p1.y, p1.z
        );
        printf(
            "(%f, %f, %f), \n",
            p2.x, p2.y, p2.z
        );
        printf(
            "(%f, %f, %f), \n",
            p3.x, p3.y, p3.z
        );
        printf(
            "(%f, %f, %f), \n",
            p4.x, p4.y, p4.z
        );
        printf(
            "(%f, %f, %f), \n",
            p5.x, p5.y, p5.z
        );
        printf(
            "(%f, %f, %f), \n",
            p6.x, p6.y, p6.z
        );
    }

    // Top-left edge point
    else if (i == 0 && j == _res1)
    {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i + 1][j];
        p3 = _prev_nodes[i + 1][j - 1];
        p4 = _prev_nodes[i + 2][j];
        p5 = _prev_nodes[i + 2][j - 1];
        p6 = _prev_nodes[i + 2][j - 2];
    }

    // Top-right edge point
    else if (i == _res1 && j == _res1)
    {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i - 1][j];
        p3 = _prev_nodes[i - 1][j - 1];
        p4 = _prev_nodes[i - 2][j];
        p5 = _prev_nodes[i - 2][j - 1];
        p6 = _prev_nodes[i - 2][j - 2];
    }

    // Bottom-right edge point
    else if (i == _res1 && j == 0)
    {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i - 1][j];
        p3 = _prev_nodes[i - 1][j + 1];
        p4 = _prev_nodes[i - 2][j];
        p5 = _prev_nodes[i - 2][j + 1];
        p6 = _prev_nodes[i - 2][j + 2];
    }

    // Bottom edge
    if (j == 0) {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i - 1][j];
        p3 = _prev_nodes[i + 1][j];
        p4 = _prev_nodes[i][j + 1];
        p5 = _prev_nodes[i - 1][j + 1];
        p6 = _prev_nodes[i + 1][j + 1];
    }

    // Right edge
    else if (i == _res1) {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i][j - 1];
        p3 = _prev_nodes[i][j + 1];
        p4 = _prev_nodes[i - 1][j];
        p5 = _prev_nodes[i - 1][j - 1];
        p6 = _prev_nodes[i - 1][j + 1];
    }

    // Top edge
    else if (j == _res1) {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i - 1][j];
        p3 = _prev_nodes[i + 1][j];
        p4 = _prev_nodes[i][j - 1];
        p5 = _prev_nodes[i - 1][j - 1];
        p6 = _prev_nodes[i + 1][j - 1];
    }

    // Left edge
    else if (i == 0) {
        p1 = _prev_nodes[i][j];
        p2 = _prev_nodes[i][j - 1];
        p3 = _prev_nodes[i][j + 1];
        p4 = _prev_nodes[i + 1][j];
        p5 = _prev_nodes[i + 1][j - 1];
        p6 = _prev_nodes[i + 1][j + 1];
    }




    // Approximate new z component
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
    const double c = 180.0 / PI;

    cramer(matrix, output, coeff);
    const double step = 0.5 / _res;


    Node new_node = p1 - vector_gradient(i, j) * step;
    
    

    double new_z = coeff[0] * new_node.x * new_node.x +
                   coeff[1] * new_node.x * new_node.y +
                   coeff[2] * new_node.y * new_node.y +
                   coeff[3] * new_node.x +
                   coeff[4] * new_node.y +
                   coeff[5];

    new_node.z = new_z;
    printf(
        "(%f, %f, %f) -> (%f, %f, %f) [%d %d]",
        p1.x, p1.y, p1.z,
        new_node.x, new_node.y, new_node.z,
        i, j
    );


    //printf("%f, ", new_z);
    double angle = atan2(new_z, step) * c;
    // printf("| %f, ", angle);

    return angle;
}

/**
 * Returns whether a point (i,j) is on any printed region.
 *
 * @brief	If the node(i,j) is on printed region.
 * @param	i	int	  The angular index of the node array.
 * @param	j	int	  The radial index of the node array.
 * @return	onreg	bool	True if the node is on the printed region.
 *                          False if the node is not on the printed region.
 */
bool Mesh::OnPrintedRegion(int i, int j)
{
    // Absolute positions
    double abs_y_val = abs(_prev_nodes[i][j].y) - half_p;

    // Region booleans
    bool r1 = (     -half_p      <= abs_y_val && abs_y_val <= 0 * w_g + 0 * w_p);
    bool r2 = (1 * w_g + 0 * w_p <= abs_y_val && abs_y_val <= 1 * w_g + 1 * w_p);
    bool r3 = (2 * w_g + 1 * w_p <= abs_y_val && abs_y_val <= 2 * w_g + 2 * w_p);
    bool r4 = (3 * w_g + 2 * w_p <= abs_y_val && abs_y_val <= 3 * w_g + 3 * w_p);

    // Return if point is on printed region
    return r1 || r2;
}



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
    int n = _res / 2;
    Node v1, v2;

    // Center
    if (i == n && j == n)
    {
        return Node(0, 0, 1);
    }

    // Top-left diagonal
    else if (i == _res1 - j && i < n)
    {  
        v1 = _prev_nodes[i - 1][j + 1] - _prev_nodes[i + 1][j - 1];  // Radius change
        v2 = _prev_nodes[i][j - 1] - _prev_nodes[i + 1][j];          // Polar angle change
    }

    // Bottom-right diagonal
    else if (i == _res1 - j && i > n)
    {  
        v1 = _prev_nodes[i + 1][j - 1] - _prev_nodes[i - 1][j + 1];  // Radius change
        v2 = _prev_nodes[i][j + 1] - _prev_nodes[i - 1][j];          // Polar angle change
    }

    // Bottom-left diagonal
    else if (i == j && i < n)
    {  
        v1 = _prev_nodes[i - 1][j - 1] - _prev_nodes[i + 1][j + 1];  // Radius change
        v2 = _prev_nodes[i + 1][j] - _prev_nodes[i][j + 1];          // Polar angle change
    }

    // Top-right diagonal
    else if (i == j && i > n)
    {
        v1 = _prev_nodes[i + 1][j + 1] - _prev_nodes[i - 1][j - 1];  // Radius change
        v2 = _prev_nodes[i - 1][j] - _prev_nodes[i][j - 1];          // Polar angle change
    }

    // Non-diagonal points
    else {
        v1 = _prev_nodes[i + 1][j] - _prev_nodes[i - 1][j];  // Normal neighbor differential
        v2 = _prev_nodes[i][j + 1] - _prev_nodes[i][j - 1];  // Normal neighbor differential
    }

    // Return normal vector
    return cross_product(v1, v2).normalize();
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
    Node v_cross = _prev_nodes[i][j - 1] + _prev_nodes[i][j + 1] + _prev_nodes[i + 1][j] + _prev_nodes[i - 1][j] - _prev_nodes[i][j] * 4;
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
    Node v1 = _prev_nodes[i - 1][j] - _prev_nodes[i][j - 1];
    Node v2 = _prev_nodes[i - 1][j] - _prev_nodes[i][j + 1];
    Node v3 = _prev_nodes[i + 1][j] - _prev_nodes[i][j - 1];
    Node v4 = _prev_nodes[i + 1][j] - _prev_nodes[i][j + 1];

    // Calculated cross products
    Node vec1 = cross_product(v1, v2);
    Node vec2 = cross_product(v3, v4);

    // Return normal pressure force
    return vector_normal(i, j) * (_pressure / 4) * (vec1.det() + vec2.det());
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





