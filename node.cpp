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

Node::Node()
{
	x = 0.0; y = 0.0; z = 0.0;
}

Node::Node(const Node& _node)
{
	x = _node.x; y = _node.y; z = _node.z;
}

Node::Node(double _x, double _y, double _z)
{
	x = _x; y = _y; z = _z;
}


Node& Node::operator= (const Node& _node)
{
	x = _node.x; y = _node.y; z = _node.z;
	return *this;
}

Node  Node::operator+ (const Node& _node)
{
	return Node(x + _node.x, y + _node.y, z + _node.z);
}

Node& Node::operator+=(const Node& _node)
{
	x += _node.x; y += _node.y; z += _node.z;
	return *this;
}

Node  Node::operator- (const Node& _node)
{
	return Node(x - _node.x, y - _node.y, z - _node.z);
}

Node& Node::operator-=(const Node& _node)
{
	x -= _node.x; y -= _node.y; z -= _node.z;
	return *this;
}

Node  Node::operator* (double _scalar)
{
	return Node(x * _scalar, y * _scalar, z * _scalar);
}

Node& Node::operator*=(double _scalar)
{
	x *= _scalar; y *= _scalar; z *= _scalar;
	return *this;
}

Node  Node::operator/ (double _scalar)
{
	assert(_scalar != 0);
	double inv_scalar = 1.0 / _scalar;
	return Node(x * inv_scalar, y * inv_scalar, z * inv_scalar);
}

Node& Node::operator/=(double _scalar)
{
	assert(_scalar != 0);
	double inv_scalar = 1.0 / _scalar;
	x *= inv_scalar; y *= inv_scalar; z *= inv_scalar;
	return *this;
}


/*
 * Normalizes the magnetude of the node.
 * @brief	Implicit node normalization.
 */
Node Node::normalize()
{
	*this = *this * f_inv_sqrt(x * x + y * y + z * z);
	return *this;
}

/*
 * Projects the node onto another node.
 * @brief	Implicit node projection.
 */
Node Node::project(Node _onto)
{
	*this = _onto * (dot_product(*this, _onto) / _onto.magnetude_squared());
	return *this;
}

/*
 * Calculates the magnetude of the node.
 * @brief	Explicit node magnetude.
 */
double Node::magnetude()
{
	return f_sqrt(x * x + y * y + z * z);
}

/*
 * Calculates the square of the magnetude of the node.
 * @brief	Explicit node square magnetude.
 */
double Node::magnetude_squared()
{
	return x * x + y * y + z * z;
}


/*
 * Calculates the dot product of 2 nodes.
 * @brief	Vector dot product.
 */
double dot_product(Node _node1, Node _node2)
{
	return _node1.x * _node2.x + _node1.y * _node2.y + _node1.z * _node2.z;
}

/*
 * Calculates the cross product of 2 nodes.
 * @brief	Vector cross product.
 */
Node cross_product(Node _node1, Node _node2)
{
	double i = _node1.y * _node2.z - _node1.z * _node2.y;
	double j = _node1.z * _node2.x - _node1.x * _node2.z;
	double k = _node1.x * _node2.y - _node1.y * _node2.x;
	return Node(i, j, k);
}