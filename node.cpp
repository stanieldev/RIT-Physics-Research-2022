//
// File: node.cpp
// Author: Stanley Goodwin
// Creation Date: 6/5/2022
// Last Modified: 7/4/2022
//
#include <assert.h>
#include <math.h>
#include "fmath.h"
#include "node.h"


// Default constructors
Node::Node()
{
	x = 0; y = 0; z = 0;
}
Node::Node(double ix, double iy, double iz)
{
	x = ix; y = iy; z = iz;
}
Node::Node(const Node& node)
{
	x = node.x; y = node.y; z = node.z;
}


// Equals
Node& Node::operator= (const Node& node)
{
	x = node.x; 
	y = node.y; 
	z = node.z;
	return *this;
}

// Addition
Node  Node::operator+ (const Node& node)
{
	return Node(x + node.x, y + node.y, z + node.z);
}
Node& Node::operator+=(const Node& node)
{
	x += node.x; 
	y += node.y; 
	z += node.z;
	return *this;
}

// Subtraction
Node  Node::operator- (const Node& node)
{
	return Node(x - node.x, y - node.y, z - node.z);
}
Node& Node::operator-=(const Node& node)
{
	x -= node.x; 
	y -= node.y; 
	z -= node.z;
	return *this;
}

// Multiplication [By Scalar]
Node  Node::operator* (double scalar)
{
	return Node(x * scalar, y * scalar, z * scalar);
}
Node& Node::operator*=(double scalar)
{
	x *= scalar;
	y *= scalar;
	z *= scalar;
	return *this;
}

// Division [By Scalar]
Node  Node::operator/ (double scalar)
{
	assert(scalar != 0);
	double inv_scalar = 1 / scalar;
	return Node(x * inv_scalar, y * inv_scalar, z * inv_scalar);
}
Node& Node::operator/=(double scalar)
{
	assert(scalar != 0);
	double inv_scalar = 1 / scalar;
	x *= inv_scalar;
	y *= inv_scalar;
	z *= inv_scalar;
	return *this;
}

// Normalization function
Node& Node::normalize() {
	double inv_magnetude = f_inv_sqrt(x * x + y * y + z * z);
	x *= inv_magnetude;
	y *= inv_magnetude;
	z *= inv_magnetude;
	return *this;
}

// Projection function
Node Node::proj(Node v)
{
	double scalar = dot_product(*this, v) / v.det2();
	return v * scalar;
}

// Characteristic functions
double Node::det()
{
	return f_sqrt(x * x + y * y + z * z);
}
double Node::det2()
{
	return x * x + y * y + z * z;
}


// Vector Arithmetic
double dot_product(Node v1, Node v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Node cross_product(Node v1, Node v2)
{
	double i = v1.y * v2.z - v1.z * v2.y;
	double j = v1.z * v2.x - v1.x * v2.z;
	double k = v1.x * v2.y - v1.y * v2.x;
	return Node(i, j, k);
}
