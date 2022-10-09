/**
 * @file	droplet.h
 * @brief	Stores the simulation's droplet parameters.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/25/2022
 * Last Modified: 7/25/2022
 */
#pragma once
#ifndef DROPLET_H
#define DROPLET_H


class Droplet {
public:
    double radius;  // Expected contact radius
    double volume;  // Expected final volume
    double vkappa;  // Ratio of volume to r^3
    double m_rad3;  // Expected contact radius cubed

    Droplet(double _contact_radius, double _volume)
    {
        radius = _contact_radius;
        volume = _volume;
        m_rad3 = radius * radius * radius;
        vkappa = volume / m_rad3;
    }
};


#endif DROPLET_H