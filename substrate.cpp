/**
 * File:	substrate.cpp
 * Brief:	Substrate-related functions.
 *
 * Author:	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 7/26/2022
 * Last Modified: 7/31/2022
 */

#include "substrate.h"

// STATUS: COMPLETE & VERIFIED
// Tests if node is on the printed region and slips
bool Substrate::slips_on_printed(Node node, double node_contact_angle)
{
    // Variable
    double abs_y_val;

    // Tests if the contact angle is large enough to not slip
    if (node_contact_angle >= m_theta_p) {
        return false;
    }
    else {
        // Declare variable values
        abs_y_val = fabs(node.y) - 0.5 * m_width_p;

        // Region booleans
        bool r1 = (       -0.5 * m_width_p       <= abs_y_val && abs_y_val <= 0 * m_width_g + 0 * m_width_p);
        bool r2 = (1 * m_width_g + 0 * m_width_p <= abs_y_val && abs_y_val <= 1 * m_width_g + 1 * m_width_p);
        /*
        bool r3 = (2 * m_width_g + 1 * m_width_p <= abs_y_val && abs_y_val <= 2 * m_width_g + 2 * m_width_p);
        bool r4 = (3 * m_width_g + 2 * m_width_p <= abs_y_val && abs_y_val <= 3 * m_width_g + 3 * m_width_p);
        bool r5 = (4 * m_width_g + 3 * m_width_p <= abs_y_val && abs_y_val <= 4 * m_width_g + 4 * m_width_p);
        bool r6 = (5 * m_width_g + 4 * m_width_p <= abs_y_val && abs_y_val <= 5 * m_width_g + 5 * m_width_p);
        bool r7 = (6 * m_width_g + 5 * m_width_p <= abs_y_val && abs_y_val <= 6 * m_width_g + 6 * m_width_p);
        bool r8 = (7 * m_width_g + 6 * m_width_p <= abs_y_val && abs_y_val <= 7 * m_width_g + 7 * m_width_p);
        */

        // Return true if point is on printed region, else false
        return r1 || r2;  // || r3 || r4 || r5 || r6 || r7 || r8
    }
}

// STATUS: COMPLETE
// Tests if node is on the non-printed region and slips
bool Substrate::slips_on_surface(Node node, double node_contact_angle)
{
    // Variable
    double abs_y_val;

    // Tests if the contact angle is large enough to not slip
    if (node_contact_angle >= m_theta_g) {
        return false;
    }
    else {
        // Declare variable values
        abs_y_val = fabs(node.y) - 0.5 * m_width_p;

        // Region booleans
        bool r1 = (       -0.5 * m_width_p       <= abs_y_val && abs_y_val <= 0 * m_width_g + 0 * m_width_p);
        bool r2 = (1 * m_width_g + 0 * m_width_p <= abs_y_val && abs_y_val <= 1 * m_width_g + 1 * m_width_p);
        /*
        bool r3 = (2 * m_width_g + 1 * m_width_p <= abs_y_val && abs_y_val <= 2 * m_width_g + 2 * m_width_p);
        bool r4 = (3 * m_width_g + 2 * m_width_p <= abs_y_val && abs_y_val <= 3 * m_width_g + 3 * m_width_p);
        bool r5 = (4 * m_width_g + 3 * m_width_p <= abs_y_val && abs_y_val <= 4 * m_width_g + 4 * m_width_p);
        bool r6 = (5 * m_width_g + 4 * m_width_p <= abs_y_val && abs_y_val <= 5 * m_width_g + 5 * m_width_p);
        bool r7 = (6 * m_width_g + 5 * m_width_p <= abs_y_val && abs_y_val <= 6 * m_width_g + 6 * m_width_p);
        bool r8 = (7 * m_width_g + 6 * m_width_p <= abs_y_val && abs_y_val <= 7 * m_width_g + 7 * m_width_p);
        */

        // Return true if point is on non-printed region, else false
        return !(r1 || r2);  // || r3 || r4 || r5 || r6 || r7 || r8
    }
}
