//
// File: node.h
// Author: Stanley Goodwin
// Creation Date: 6/5/2022
// Last Modified: 6/20/2022
//
#pragma once
#ifndef NODE_H
#define NODE_H


/*
The class that handles the node position vectors and operations involving 3-vectors.
*/
class Node {
public:
	// Stored values
	double x;
	double y;
	double z;

	// Default constuctors
	Node();
	Node(double ix, double iy, double iz);
	Node(const Node& node);

	// C++ arithmetic operations
	Node& operator= (const Node& node);
	Node  operator+ (const Node& node);
	Node& operator+=(const Node& node);
	Node  operator- (const Node& node);
	Node& operator-=(const Node& node);
	Node  operator* (double scalar);
	Node& operator*=(double scalar);
	Node  operator/ (double scalar);
	Node& operator/=(double scalar);

	// Other vector operations
	Node& normalize();
	double det();
	double det2();
	Node proj(Node v);
};


// Node functions
double dot_product(Node v1, Node v2);
Node cross_product(Node v1, Node v2);



#endif NODE_H