/**
 * @file	droplet.h
 * @brief	Stores droplet characteristics in a single class.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/16/2022
 * Last Modified: 7/5/2022
 */
#pragma once
#ifndef DROPLET_H
#define DROPLET_H


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
    double radius   = 1.9407025E-3;
    double volume   = 3.0E-9;
    double pressure = 1.0E+5;

    double κ = volume / (radius * radius * radius);
};


#endif DROPLET_H