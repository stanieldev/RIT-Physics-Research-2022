/*
 * File:	droplet.cpp
 * Author:	Stanley Goodwin
 * Stores the simulation's expected droplet parameters.
 */
#include "droplet.hpp"

/*
 * Definition of constructing a droplet-type object.
 * @brief	Droplet contructor definition.
 * @param	_contact_radius double	The contact radius of the droplet.
 * @param	_volume         double	The volume of the droplet.
 */
Droplet::Droplet(double _contact_radius, double _volume)
{
    contact_radius = _contact_radius;
    contact_radius_cubed = _contact_radius * _contact_radius * _contact_radius;
    volume = _volume;
    volume_ratio = volume / contact_radius_cubed;
}

/*
 * Rescales the volume ratio by a scale factor.
 * @brief   Rescales droplet volume.
 * @param	_volume_factor  double	The factor to scale by.
 */
void Droplet::rescale(double _volume_factor)
{
    volume_ratio *= _volume_factor;
}