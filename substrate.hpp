/*
 * File:	substrate.hpp
 * Author:	Stanley Goodwin
 * Stores the simulation's substrate parameters.
 */
#pragma once
#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "node.hpp"
#include <math.h>

struct Substrate {
    double printed_receding_angle, printed_receding_width;
    double unprinted_receding_angle, unprinted_receding_width;
    double kp, kg;

    Substrate(
        double _printed_receding_angle,
        double _printed_receding_width,
        double _unprinted_receding_angle,
        double _unprinted_receding_width
    );

    bool slips_on_printed(Node node, double contact_angle);
    bool slips_on_unprinted(Node node, double contact_angle);
};

#endif SUBSTRATE_H