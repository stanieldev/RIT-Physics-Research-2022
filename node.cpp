/*
 * File:	node.cpp
 * Author:	Stanley Goodwin
 *
 * Stores all the necessary data of the 3-component
 * vector, as well as the operations between node-node
 * and node-scalar. Includes some other vector functions
 * like projection, normalization, and determinant.
 */
#include <assert.h>
#include "fmath.hpp"
#include "node.hpp"

/*
 * Initializes a node to (0, 0, 0) if no input specified.
 * @brief	Default constructor.
 */
Node::Node()
{
	x = 0; y = 0; z = 0;
}

/*
 * Initializes a node using user coordinates.
 * @brief	Coordinate-specified constructor.
 * @param	ix	double	The initial x coordinate.
 * @param	iy	double	The initial y coordinate.
 * @param	iz	double	The initial z coordinate.
 */
Node::Node(double ix, double iy, double iz)
{
	x = ix; y = iy; z = iz;
}

/*
 * Initializes a node using another node as reference.
 * @brief	Node-specified constructor.
 * @param	node	Node	The node to be copies.
 */
Node::Node(const Node& node)
{
	x = node.x; y = node.y; z = node.z;
}


/*
 * Copies the contents of a node to another node.
 * @brief	Node equality.
 * @param	node	Node	The node to be copied.
 * @return	result	Node	A new copy of the node.
 */
Node& Node::operator= (const Node& node)
{
	x = node.x; y = node.y; z = node.z;
	return *this;
}

/*
 * Adds 2 node's components together and returns the result as another node.
 * @brief	Node addition.
 * @param	node	Node	The node to be added.
 * @return	result	Node	The sum of the nodes.
 */
Node  Node::operator+ (const Node& node)
{
	return Node(x + node.x, y + node.y, z + node.z);
}

/*
 * Adds the input node's components to the original node.
 * @brief	Implicit node addition.
 * @param	node	Node	The node to be added to the other.
 */
Node& Node::operator+=(const Node& node)
{
	x += node.x; y += node.y; z += node.z;
	return *this;
}

/*
 * Subtracts 2 node's components together and returns the result as another node.
 * @brief	Node subtraction.
 * @param	node	Node	The node to be subtracted.
 * @return	result	Node	The difference of the nodes.
 */
Node  Node::operator- (const Node& node)
{
	return Node(x - node.x, y - node.y, z - node.z);
}

/*
 * Subtracts the input node's components from the original node.
 * @brief	Implicit node subtraction.
 * @param	node	Node	The node to be subtracted to the other.
 */
Node& Node::operator-=(const Node& node)
{
	x -= node.x; y -= node.y; z -= node.z;
	return *this;
}

/*
 * Multiplies the node's components by a scalar and returns the result as another node.
 * @brief	Scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 * @return	result	Node	The scaled node.
 */
Node  Node::operator* (double scalar)
{
	return Node(x * scalar, y * scalar, z * scalar);
}

/*
 * Multiplies the node's components by a scalar.
 * @brief	Implicit scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 */
Node& Node::operator*=(double scalar)
{
	x *= scalar; y *= scalar; z *= scalar;
	return *this;
}

/*
 * Divides the node's components by a scalar and returns the result as another node.
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

/*
 * Divides the node's components by a scalar.
 * @brief	Implicit scalar multiplication.
 * @param	scalar	double	The scalar to multiple the node by.
 */
Node& Node::operator/=(double scalar)
{
	assert(scalar != 0);
	double inv_scalar = 1 / scalar;
	x *= inv_scalar; y *= inv_scalar; z *= inv_scalar;
	return *this;
}


/*
 * Normalizes the magnetude of the node.
 * @brief	Node normalization.
 */
Node& Node::normalize() {
	double inv_magnetude = f_inv_sqrt(x * x + y * y + z * z);
	x *= inv_magnetude; y *= inv_magnetude; z *= inv_magnetude;
	return *this;
}

/*
 * Projects the original node onto another node.
 * @brief	Returns a new projected node.
 * @param	onto	Node	The node to do the projection onto.
 */
void Node::project(Node onto)
{
	*this = onto * (dot_product(*this, onto) / onto.magnetude_squared());
}

/*
 * Calculates the magnetude of a node.
 * @brief	Node determinant / magnetude.
 * @return	magn	double	The magnetude of the node.
 */
double Node::magnetude()
{
	return f_sqrt(x * x + y * y + z * z);
}

/*
 * Calculates the square of the magnetude of a node.
 * @brief	Node determinant / magnetude squared.
 * @return	magn2	double	The square magnetude of the node.
 */
double Node::magnetude_squared()
{
	return x * x + y * y + z * z;
}


/*
 * Calculates the dot product of 2 nodes.
 * @brief	Vector dot product.
 * @param	node1	Node	The first node.
 * @param	node2	Node	The second node.
 * @return	dotp	double	The vector dot product.
 */
double dot_product(Node node1, Node node2)
{
	return node1.x * node2.x + node1.y * node2.y + node1.z * node2.z;
}

/*
 * Calculates the cross product of 2 nodes.
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

/*
 * Projects the original node onto an input node and returns the projected vector.
 * @brief	Returns a new projected node.
 * @param	node	Node	The node to project.
 * @param	onto	Node	The node to do the projection onto.
 * @return	projected_node	Node	The vector projection.
 */
Node project(Node node, Node onto)
{
	return onto * (dot_product(node, onto) / onto.magnetude_squared());
}