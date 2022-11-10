// FINISHED
/*
 * File:	node.cpp
 * Author:	Stanley Goodwin
 *
 * Stores all the necessary data of the 3-component vector,
 * Allows operations between node-node and node-scalar.
 * Includes some other vector functions like projection, normalization, and magnetude.
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
 * Initializes a node using input coordinates.
 * @brief	Coordinate-specified constructor.
 * @param	_x	double	The initial x-coordinate.
 * @param	_y	double	The initial y-coordinate.
 * @param	_z	double	The initial z-coordinate.
 */
Node::Node(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}

/*
 * Initializes a node using another node as reference.
 * @brief	Node-specified constructor.
 * @param	_node	Node	The copy of the input node.
 */
Node::Node(const Node& _node)
{
	x = _node.x; y = _node.y; z = _node.z;
}


/*
 * Node component-wise assignment.
 * @brief	Node assignment.
 * @param	_node	Node
 */
Node& Node::operator= (const Node& _node)
{
	x = _node.x; y = _node.y; z = _node.z;
	return *this;
}

/*
 * Node component-wise addition.
 * @brief	Explicit node addition.
 * @param	_node	Node
 */
Node  Node::operator+ (const Node& _node)
{
	return Node(x + _node.x, y + _node.y, z + _node.z);
}

/*
 * Node component-wise implicit addition.
 * @brief	Implicit node addition.
 * @param	_node	Node
 */
Node& Node::operator+=(const Node& _node)
{
	x += _node.x; y += _node.y; z += _node.z;
	return *this;
}

/*
 * Node component-wise subtraction.
 * @brief	Explicit node subtraction.
 * @param	_node	Node
 */
Node  Node::operator- (const Node& _node)
{
	return Node(x - _node.x, y - _node.y, z - _node.z);
}

/*
 * Node component-wise implicit subtraction.
 * @brief	Implicit node subtraction.
 * @param	_node	Node
 */
Node& Node::operator-=(const Node& _node)
{
	x -= _node.x; y -= _node.y; z -= _node.z;
	return *this;
}

/*
 * Node component-wise scalar multiplication.
 * @brief	Explicit node scalar multiplication.
 * @param	_scalar	double
 */
Node  Node::operator* (double _scalar)
{
	return Node(x * _scalar, y * _scalar, z * _scalar);
}

/*
 * Node component-wise implicit scalar multiplication.
 * @brief	Implicit node scalar multiplication.
 * @param	_scalar	double
 */
Node& Node::operator*=(double _scalar)
{
	x *= _scalar; y *= _scalar; z *= _scalar;
	return *this;
}

/*
 * Node component-wise scalar division.
 * @brief	Explicit node scalar division.
 * @param	_scalar	double
 */
Node  Node::operator/ (double _scalar)
{
	assert(_scalar != 0);
	double inv_scalar = 1 / _scalar;
	return Node(x * inv_scalar, y * inv_scalar, z * inv_scalar);
}

/*
 * Node component-wise implicit scalar division.
 * @brief	Implicit node scalar division.
 * @param	_scalar	double
 */
Node& Node::operator/=(double _scalar)
{
	assert(_scalar != 0);
	double inv_scalar = 1 / _scalar;
	x *= inv_scalar; y *= inv_scalar; z *= inv_scalar;
	return *this;
}


/*
 * Normalizes the magnetude of the node.
 * @brief	Implicit node normalization.
 */
void Node::normalize()
{
	*this = *this * f_inv_sqrt(x * x + y * y + z * z);
}

/*
 * Projects the node onto another node.
 * @brief	Implicit node projection.
 * @param	_onto	Node
 */
void Node::project(Node _onto)
{
	*this = _onto * (dot_product(*this, _onto) / _onto.magnetude_squared());
}

/*
 * Calculates the magnetude of the node.
 * @brief	Explicit node magnetude.
 * @return	magn	double
 */
double Node::magnetude()
{
	return f_sqrt(x * x + y * y + z * z);
}

/*
 * Calculates the square of the magnetude of the node.
 * @brief	Explicit node square magnetude.
 * @return	magn2	double
 */
double Node::magnetude_squared()
{
	return x * x + y * y + z * z;
}


/*
 * Normalizes the magnetude of the node.
 * @brief	Explicit node normalization.
 */
Node normalize(Node _node)
{
	return _node * f_inv_sqrt(_node.x * _node.x + _node.y * _node.y + _node.z * _node.z);
}

/*
 * Projects the node onto another node.
 * @brief	Explicit node projection.
 * @param	_onto	Node
 */
Node project(Node _node, Node _onto)
{
	return _onto * (dot_product(_node, _onto) / _onto.magnetude_squared());
}

/*
 * Calculates the dot product of 2 nodes.
 * @brief	Vector dot product.
 * @param	_node1	Node	The first node.
 * @param	_node2	Node	The second node.
 * @return	dotp	double
 */
double dot_product(Node _node1, Node _node2)
{
	return _node1.x * _node2.x + _node1.y * _node2.y + _node1.z * _node2.z;
}

/*
 * Calculates the cross product of 2 nodes.
 * @brief	Vector cross product.
 * @param	_node1	Node	The first node.
 * @param	_node2	Node	The second node.
 * @return	crossp	Node
 */
Node cross_product(Node _node1, Node _node2)
{
	double i = _node1.y * _node2.z - _node1.z * _node2.y;
	double j = _node1.z * _node2.x - _node1.x * _node2.z;
	double k = _node1.x * _node2.y - _node1.y * _node2.x;
	return Node(i, j, k);
}