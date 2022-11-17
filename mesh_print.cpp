/**
 * @file	mesh_print.cpp
 * @brief	Mesh printing functions.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/18/2022
 * Last Modified: 7/31/2022
 */
//#include <fstream>
//#include <string>
//#include "mesh.hpp"
//
//
///**
// * Prints the contents of the node arrays to a text file.
// * @brief   Save mesh's nodes to file.
// */
//void Mesh::fprint_nodes()
//{
//    // Variables
//    std::ofstream file_x, file_y, file_z;
//    std::string output_folder = "C:\\Users\\sfg99\\Code\\Summer Research\\Matlab\\generated_data";
//
//    std::string file_name = "nodes";
//    std::string file_end = ".txt";
//    file_name += "_" + std::to_string(iterations_run);
//    if (current_iter != 0) { file_name += "_" + std::to_string(current_iter); }
//
//
//    // Open file streams
//    file_x.open(output_folder + "\\" + file_name + "_x" + file_end);
//    file_y.open(output_folder + "\\" + file_name + "_y" + file_end);
//    file_z.open(output_folder + "\\" + file_name + "_z" + file_end);
//
//    // Print all nodes to separate files
//    Node _node;
//    for (int i = 0; i < m_res; i++)
//    {
//        for (int j = 0; j < m_res; j++)
//        {
//            // Set current node
//            _node = m_current_nodes[i][j];  //  * _droplet.contact_radius
//
//            // Send data to the files
//            file_x << _node.x << " ";
//            file_y << _node.y << " ";
//            file_z << _node.z << " ";
//        }
//
//        // Add new line
//        file_x << "\n";
//        file_y << "\n";
//        file_z << "\n";
//    }
//
//    // Close file streams
//    file_x.close();
//    file_y.close();
//    file_z.close();
//}








//void Mesh::fprint_volume()
//{
//
//};
//void Mesh::fprint_pressure()
//{
//
//};
//void Mesh::fprint_gamma()
//{
//
//};
//
//void Mesh::cprint_nodes()
//{
//
//};
//void Mesh::cprint_volume()
//{
//
//};
//void Mesh::cprint_pressure()
//{
//    
//};
//void Mesh::cprint_gamma()
//{
//
//};