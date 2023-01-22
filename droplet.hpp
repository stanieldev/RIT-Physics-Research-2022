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

    Droplet();
    Droplet(const Droplet& _droplet);
    Droplet(double _contact_radius, double _volume);
};

#endif DROPLET_H