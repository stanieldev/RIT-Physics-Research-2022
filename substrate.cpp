/*
 * File:	substrate.cpp
 * Author:	Stanley Goodwin
 * Stores the simulation's substrate characteristics.
 */
#include <math.h>
#include "substrate.hpp"

Substrate::Substrate()
{
    printed_receding_angle = 0.0;
    printed_region_width = 0.0;
    printed_invtan = 0.0;
    unprinted_receding_angle = 0.0;
    unprinted_region_width = 0.0;
    unprinted_invtan = 0.0;
}

Substrate::Substrate(const Substrate& _substrate)
{
    printed_receding_angle = _substrate.printed_receding_angle;
    printed_region_width   = _substrate.printed_region_width;
    printed_invtan         = _substrate.printed_invtan;
    unprinted_receding_angle = _substrate.unprinted_receding_angle;
    unprinted_region_width   = _substrate.unprinted_region_width;
    unprinted_invtan         = _substrate.unprinted_invtan;
}

Substrate::Substrate(
    double _printed_receding_angle_degrees,
    double _printed_region_width,
    double _unprinted_receding_angle_degrees,
    double _unprinted_region_width
) {
    constexpr auto DEG_TO_RAD = 3.141593265358979323846 / 180.0;
    printed_receding_angle = _printed_receding_angle_degrees * DEG_TO_RAD;
    printed_region_width = _printed_region_width;
    printed_invtan = 1.0 / tan(printed_receding_angle);
    unprinted_receding_angle = _unprinted_receding_angle_degrees * DEG_TO_RAD;
    unprinted_region_width = _printed_region_width;
    unprinted_invtan = 1.0 / tan(unprinted_receding_angle);
}


/*
 * Tests if node is on the printed region and slips.
 * @brief	Node slip function.
 */
bool Substrate::slips_on_printed(Node _node, double _contact_angle)
{
    // Tests if the contact angle is large enough to not slip
    if (_contact_angle >= printed_receding_angle)
        return false;

    // Declare variable values
    double abs_y_val = fabs(_node.y) - 0.5 * printed_region_width;

    // Region booleans
    bool r1 = (                          -0.5 * printed_region_width <= abs_y_val && abs_y_val <= 0 * unprinted_region_width + 0 * printed_region_width);
    bool r2 = (1 * unprinted_region_width + 0 * printed_region_width <= abs_y_val && abs_y_val <= 1 * unprinted_region_width + 1 * printed_region_width);

    // Return true if point is on printed region, else false
    return r1 || r2;
}

/*
 * Tests if node is on the unprinted region and slips.
 * @brief	Node slip function.
 */
bool Substrate::slips_on_unprinted(Node _node, double _contact_angle)
{
    // Tests if the contact angle is large enough to not slip
    if (_contact_angle >= unprinted_receding_angle)
        return false;
    
    // Declare variable values
    double abs_y_val = fabs(_node.y) - 0.5 * printed_region_width;

    // Region booleans
    bool r1 = (                          -0.5 * printed_region_width <= abs_y_val && abs_y_val <= 0 * unprinted_region_width + 0 * printed_region_width);
    bool r2 = (1 * unprinted_region_width + 0 * printed_region_width <= abs_y_val && abs_y_val <= 1 * unprinted_region_width + 1 * printed_region_width);

    // Return true if point is on non-printed region, else false
    return !(r1 || r2);
}