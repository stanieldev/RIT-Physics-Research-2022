//
// File: main.cpp
// Author: Stanley Goodwin
// Creation Date: 9/17/2021
// Last Modified: 7/4/2022
// Credit to Kara Maki for skeleton code.
//
#include <iostream>
#include "droplet.h"
#include "fmath.h"
#include "mesh.h"


int main(int argc, const char* argv[])
{
    // Create simulation objects
    Droplet droplet;
    Mesh mesh(droplet, 121, 1200);

    // Mesh functions
    mesh.InitializeNodes();
    mesh.Iterate();

    return 0;
}
