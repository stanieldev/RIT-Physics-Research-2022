/**
 * @file	mesh_print.cpp
 * @brief	Mesh printing functions.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/18/2022
 * Last Modified: 7/18/2022
 */
#include <fstream>
#include <string>
#include "mesh.h"



 /*
 * Function description:
 *   Prints the contents of the node arrays to a text file.
 *
 * Notes:
 *   Currently only works with an odd number of points.
 *   New contact angle stuff (TODO)
 */
void Mesh::fprint_nodes()
{
    // File variables
    std::ofstream fileX, fileY, fileZ;
    std::string directory = "C:\\Users\\sfg99\\Code\\Summer Research\\Matlab\\";//
    std::string file_name = "Xitwocompleted41";
    std::string file_end = ".txt";

    // Open file streams
    fileX.open(directory + file_name + "X" + std::to_string(m_iteration) + file_end);
    fileY.open(directory + file_name + "Y" + std::to_string(m_iteration) + file_end);
    fileZ.open(directory + file_name + "Z" + std::to_string(m_iteration) + file_end);

    // Print all nodes to separate files
    Node _node;
    for (int i = 0; i < m_res; i++)
    {
        for (int j = 0; j < m_res; j++)
        {
            // Set current node
            _node = m_curr_nodes[i][j];  //  * _droplet.contact_radius

            // Send data to the files
            fileX << _node.x << " ";
            fileY << _node.y << " ";
            fileZ << _node.z << " ";
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


void Mesh::fprint_volume()
{

};
void Mesh::fprint_pressure()
{

};
void Mesh::fprint_gamma()
{

};

void Mesh::cprint_nodes()
{

};
void Mesh::cprint_volume()
{

};
void Mesh::cprint_pressure()
{
    
};
void Mesh::cprint_gamma()
{

};