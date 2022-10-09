/**
 * File:	node.cpp
 * Brief:	Node-related functions.
 *
 * Author:	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/5/2022
 * Last Modified: 7/31/2022
 */
#include <assert.h>
#include "fmath.h"
#include "node.h"


 // Status: COMPLETE & VERIFIED
/**
 * Initializes a node to (0, 0, 0) if no input specified.
 * @brief	Default constructor.
 */
Node::Node()
{
	x = 0;
	y = 0;
	z = 0;
}

// Status: COMPLETE & VERIFIED
/**
 * Initializes a node using user coordinates.
 * 
 * @brief	Coordinate-specified constructor.
 * @param	ix	double	The initial x coordinate.
 * @param	iy	double	The initial y coordinate.
 * @param	iz	double	The initial z coordinate.
 */
Node::Node(double ix, double iy, double iz)
{
	x = ix;
	y = iy;
	z = iz;
}

// Status: COMPLETE & VERIFIED
/**
 * Initializes a node using another node as reference.
 *
 * @brief	Node-specified constructor.
 * @param	node	Node	The node to be copies.
 */
Node::Node(const Node& node)
{
	x = node.x;
	y = node.y;
	z = node.z;
}

// Status: COMPLETE & VERIFIED
/**
 * Copies the contents of a node to another node.
 *
 * @brief	Node equality.
 * @param	node	Node	The node to be copied.
 * @return	result	Node	A new copy of the node.
 */
Node& Node::operator= (const Node& node)
{
	x = node.x; 
	y = node.y; 
	z = node.z;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Adds 2 node's components together and 
 * returns the result as another node.
 *
 * @brief	Node addition.
 * @param	node	Node	The node to be added.
 * @return	result	Node	The sum of the nodes.
 */
Node  Node::operator+ (const Node& node)
{
	return Node(x + node.x, y + node.y, z + node.z);
}

// Status: COMPLETE & VERIFIED
/**
 * Adds the input node's components to the original node.
 *
 * @brief	Implicit node addition.
 * @param	node	Node	The node to be added to the other.
 */
Node& Node::operator+=(const Node& node)
{
	x += node.x; 
	y += node.y; 
	z += node.z;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Subtracts 2 node's components together and
 * returns the result as another node.
 * @brief	Node subtraction.
 * @param	node	Node	The node to be subtracted.
 * @return	result	Node	The difference of the nodes.
 */
Node  Node::operator- (const Node& node)
{
	return Node(x - node.x, y - node.y, z - node.z);
}

// Status: COMPLETE & VERIFIED
/**
 * Subtracts the input node's components from the original node.
 *
 * @brief	Implicit node subtraction.
 * @param	node	Node	The node to be subtracted to the other.
 */
Node& Node::operator-=(const Node& node)
{
	x -= node.x; 
	y -= node.y; 
	z -= node.z;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Multiplies the node's components by a scalar 
 * and returns the result as another node.
 *
 * @brief	Scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 * @return	result	Node	The scaled node.
 */
Node  Node::operator* (double scalar)
{
	return Node(x * scalar, y * scalar, z * scalar);
}

// Status: COMPLETE & VERIFIED
/**
 * Multiplies the node's components by a scalar.
 *
 * @brief	Implicit scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 */
Node& Node::operator*=(double scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Divides the node's components by a scalar
 * and returns the result as another node.
 *
 * @brief	Scalar division.
 * @param	scalar	double	The scalar to divide the node by.
 * @return	result	Node	The scaled node.
 */
Node  Node::operator/ (double scalar)
{
	assert(scalar != 0);
	double inv_scalar = 1 / scalar;
	return Node(x * inv_scalar, y * inv_scalar, z * inv_scalar);
}

// Status: COMPLETE & VERIFIED
/**
 * Divides the node's components by a scalar.
 *
 * @brief	Implicit scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 */
Node& Node::operator/=(double scalar)
{
	assert(scalar != 0);
	double inv_scalar = 1 / scalar;
	x *= inv_scalar;
	y *= inv_scalar;
	z *= inv_scalar;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Normalizes the magnetude of the node.
 *
 * @brief	Node normalization.
 */
Node& Node::normalize() {
	double inv_magnetude = f_inv_sqrt(x * x + y * y + z * z);
	x *= inv_magnetude;
	y *= inv_magnetude;
	z *= inv_magnetude;
	return *this;
}

// Status: COMPLETE & VERIFIED
/**
 * Projects the original node onto an input node
 * and returns the projected vector.
 *
 * @brief	Node projection.
 * @param	node	Node	The node to do the projection onto.
 * @return	projec	Node	The vector projection.
 */
Node Node::proj(Node node)
{
	double scalar = dot_product(*this, node) / node.det2();
	return node * scalar;
}

// Status: COMPLETE & VERIFIED
/**
 * Calculates the magnetude of a node.
 *
 * @brief	Node determinant / magnetude.
 * @return	magn	double	The magnetude of the node.
 */
double Node::det()
{
	return f_sqrt(x * x + y * y + z * z);
}

// Status: COMPLETE & VERIFIED
/**
 * Calculates the square of the magnetude of a node.
 *
 * @brief	Node determinant / magnetude squared.
 * @return	magn2	double	The square magnetude of the node.
 */
double Node::det2()
{
	return x * x + y * y + z * z;
}

// Status: COMPLETE & VERIFIED
/**
 * Calculates the dot product of 2 nodes.
 * 
 * @brief	Vector dot product.
 * @param	node1	Node	The first node.
 * @param	node2	Node	The second node.
 * @return	dotp	double	The vector dot product.
 */
double dot_product(Node node1, Node node2)
{
	return node1.x * node2.x + node1.y * node2.y + node1.z * node2.z;
}

// Status: COMPLETE & VERIFIED
/**
 * Calculates the cross product of 2 nodes.
 *
 * @brief	Vector cross product.
 * @param	node1	Node	The first node.
 * @param	node2	Node	The second node.
 * @return	crossp	Node	The vector cross product.
 */
Node cross_product(Node node1, Node node2)
{
	double i = node1.y * node2.z - node1.z * node2.y;
	double j = node1.z * node2.x - node1.x * node2.z;
	double k = node1.x * node2.y - node1.y * node2.x;
	return Node(i, j, k);
}
