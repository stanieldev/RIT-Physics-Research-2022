//
// File: main.cpp
// Author: Stanley Goodwin
// Creation Date: 9/17/2021
// Last Modified: 6/23/2022
// Credit to Kara Maki for skeleton code.
//
#include <iostream>
#include "mesh.h"
#include "droplet.h"
#include "fmath.h"


int main(int argc, const char* argv[])
{
    // Create objects
    Mesh mesh;
    Droplet droplet;

    // Mesh initialization
    mesh.InitializeNodes();

    // Run mesh iterations
    mesh.Iterate();

    // Print results to files
    mesh.vprintCurrentMassInformation("C:\\Users\\sfg99\\Code\\Summer Research\\Matlab\\data_generated\\Xitwocompleted41", droplet);

    return 0;
}
