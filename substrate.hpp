/*
 * File:	substrate.hpp
 * Author:	Stanley Goodwin
 * Stores the simulation's substrate characteristics.
 */
#pragma once
#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "node.hpp"

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

    bool slips_on_printed(Node _node, double _contact_angle);
    bool slips_on_unprinted(Node _node, double _contact_angle);
};

#endif SUBSTRATE_H