/**
 * @file	substrate.h
 * @brief	Stores the simulation's substrate parameters.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/25/2022
 * Last Modified: 7/25/2022
 */
#pragma once
#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "node.h"


class Substrate {
private:
    double m_theta_p;  // Receding contact angle of the printed regions
    double m_theta_g;  // Receding contact angle of the gap regions
    double m_width_p;  // Width of the printed regions
    double m_width_g;  // Width of the gap regions

    // TEMPORARY
public:
    double k1 = 1.0 / tan(m_theta_p);
    double k2 = 1.0 / tan(m_theta_g);

public:
    bool slips_on_printed(Node node, double node_contact_angle);  // If node slips on printed region
    bool slips_on_surface(Node node, double node_contact_angle);  // If node slips on non-printed region

    Substrate(
        double _theta_p, double _width_p,
        double _theta_g, double _width_g
    ) {
        m_theta_p = _theta_p;
        m_width_p = _width_p;
        m_theta_g = _theta_g;
        m_width_g = _width_g;
    }
};


#endif SUBSTRATE_H