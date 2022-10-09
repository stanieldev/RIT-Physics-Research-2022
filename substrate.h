/**
 * File:	substrate.h
 * Brief:	Stores the simulation's substrate parameters.
 *
 * Author:	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/26/2022
 * Last Modified: 7/31/2022
 */
#pragma once
#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "node.h"
#include <math.h>


class Substrate {
private:
    double m_theta_p;  // Receding contact angle of the printed regions
    double m_theta_g;  // Receding contact angle of the gap regions
    double m_width_p;  // Width of the printed regions
    double m_width_g;  // Width of the gap regions

public:
    double kp;
    double kg;
    bool slips_on_printed(Node node, double node_contact_angle);  // If node slips on printed regions
    bool slips_on_surface(Node node, double node_contact_angle);  // If node slips on gap regions

    Substrate()  // Not Needed
    {
        m_theta_p = 0.0;
        m_width_p = 0.0;
        m_theta_g = 0.0;
        m_width_g = 0.0;

        kp = 1.0;
        kg = 1.0;
    }
    Substrate(
        double _theta_p, double _width_p,
        double _theta_g, double _width_g
    ) {
        m_theta_p = _theta_p;
        m_width_p = _width_p;
        m_theta_g = _theta_g;
        m_width_g = _width_g;

        kp = 1.0 / tan(m_theta_p);
        kg = 1.0 / tan(m_theta_g);
    }
};


#endif SUBSTRATE_H