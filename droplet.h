//
// File: droplet.h
// Author: Stanley Goodwin
// Creation Date: 6/16/2022
// Last Modified: 7/4/2022
//
#pragma once
#ifndef DROPLET_H
#define DROPLET_H


/*
The class that handles the droplet's characteristics.
*/
class Droplet {
public:
    // Stored values (change for different droplets)
    double radius = 1.9407025E-3;
    double volume = 3.0E-9;
    double pressure = 1.0E5;

    // Calculated values
    double κ = volume / (radius * radius * radius);
};


#endif DROPLET_H