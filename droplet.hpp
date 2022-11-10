// FINISHED
/*
 * File:	droplet.hpp
 * Author:	Stanley Goodwin
 * Stores the simulation's droplet parameters.
 */
#pragma once
#ifndef DROPLET_H
#define DROPLET_H

struct Droplet {
    double contact_radius;
    double contact_radius_cubed;
    double volume;
    double volume_ratio;

    Droplet(double _contact_radius, double _volume);
    void rescale(double _volume_factor);
};

#endif DROPLET_H