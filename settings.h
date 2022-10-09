/**
 * @file	settings.h
 * @brief	Stores all modifiable user parameters.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/5/2022
 * Last Modified: 7/6/2022
 */
#pragma once
#ifndef SETTINGS_H
#define SETTINGS_H

#define PI 3.141593265358979323846


/**
 * Droplet characteristics class.
 *
 * Stores all the necessary data of the droplet
 * that the fishnet algorithm is attempting to
 * approximate to match the data within this class.
 *
 */
class Droplet {
public:
    double radius = 0.183469627736;  //1.9407025E-3;
    double volume = 0.0013492135069; //3.0E-9;
    //double pressure = 1.0E+5;

    double κ = volume / (radius * radius * radius);
};


/**
 * Simulation parameters class.
 *
 * The control panel for the simulation
 * variables and constants.
 *
 */
class MeshSettings {
protected:
    // Simulation settings
    int _nλ = 12000;  // Total mesh iterations per frame
    int _res = 121;  // The number of angular & radial subdivisions

    // Droplet settings
    Droplet _droplet;

    // Substrate settings
    double w_g = 0.5 * 0.001 / _droplet.radius;
    double w_p = 0.5 * 0.001 / _droplet.radius;
    double _θi = 78 * PI / 180;  // Initial contact angle
    double θ_c = -77 * PI / 180;  // Printed region contact angle

    // Mesh strength parameters
    double µ = 8.0e2;    // Frictional decay constant
    double σ = 6.0e1;    // Coefficient of surface tension
    double τ = 1.0e2;    // Tangential "spring" constant
    double δ = 1.0e4;    // relaxation factor for the pressure
    double α = 0.5;      // Boundary relaxation factor (coding alpha)
    double β = 1.0 - α;  // Conjugate relaxation factor (coding beta)
};


#endif SETTINGS_H