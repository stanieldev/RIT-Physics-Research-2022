/**
 * @file	settings.h
 * @brief	Stores all modifiable user parameters.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/5/2022
 * Last Modified: 7/20/2022
 */
#pragma once
#ifndef SETTINGS_H
#define SETTINGS_H

#define PI 3.141593265358979323846
#define experimental_droplet false




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
#if experimental_droplet
    const double drop_radius = 9.781e-2;  // Contact radius of the droplet
    const double drop_volume = 1.416e-3;  // Expected volume of the droplet
#endif



    // Substrate settings [To be added to on_printed_region()]
#if experimental_droplet
    double w_g = 0.5 * 0.01 / 9.781e-2;
    double w_p = 0.5 * 0.01 / 9.781e-2;
#else experimental_droplet
    double w_g = 0.5 * 0.001 / 1.9407025E-3;
    double w_p = 0.5 * 0.001 / 1.9407025E-3;
#endif





    // Mesh strength parameters
    const double µ = 8.0e2;    // Frictional decay constant
    const double σ = 6.0e1;    // Coefficient of surface tension
    const double τ = 1.0e2;    // Tangential "spring" constant
          double δ = 1.0e4;    // relaxation factor for the pressure
    const double α = 0.5;      // Boundary relaxation factor (coding alpha)
    const double β = 1.0 - α;  // Conjugate relaxation factor (coding beta)
};


#endif SETTINGS_H