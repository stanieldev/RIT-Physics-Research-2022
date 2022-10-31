/*
 * File:	droplet.h
 * Author:	Stanley Goodwin
 * Stores the simulation's droplet parameters.
 */
#pragma once
#ifndef DROPLET_H
#define DROPLET_H

struct Droplet {
    double radius;   // Expected contact radius
    double radius3;  // Expected contact radius cubed
    double volume;   // Expected final volume
    double vkappa;   // Ratio of volume to r^3

    Droplet(double _contact_radius, double _volume)
    {
        radius = _contact_radius;
        radius3 = radius * radius * radius;
        volume = _volume;
        vkappa = volume / radius3;
    }
};

#endif DROPLET_H