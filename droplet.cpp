/*
 * File:	droplet.cpp
 * Author:	Stanley Goodwin
 * Stores the simulation's droplet parameters.
 */
#include "droplet.hpp"

Droplet::Droplet()
{
    contact_radius = 0.0;
    contact_radius_cubed = 0.0;
    volume = 0.0;
    volume_ratio = 1.0;
}

Droplet::Droplet(const Droplet& _droplet)
{
    contact_radius = _droplet.contact_radius;
    contact_radius_cubed = _droplet.contact_radius_cubed;
    volume = _droplet.volume;
    volume_ratio = _droplet.volume_ratio;
}

Droplet::Droplet(double _contact_radius, double _volume)
{
    contact_radius = _contact_radius;
    contact_radius_cubed = _contact_radius * _contact_radius * _contact_radius;
    volume = _volume;
    volume_ratio = _volume / contact_radius_cubed;
}