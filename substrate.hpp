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
    double printed_receding_angle, printed_region_width, printed_invtan;
    double unprinted_receding_angle, unprinted_region_width, unprinted_invtan;

    Substrate();
    Substrate(const Substrate& _substrate);
    Substrate(
        double _printed_receding_angle_degrees,
        double _printed_region_width,
        double _unprinted_receding_angle_degrees,
        double _unprinted_region_width
    );

    bool slips_on_printed(Node _node, double _contact_angle);
    bool slips_on_unprinted(Node _node, double _contact_angle);
};

#endif SUBSTRATE_H