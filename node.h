//
// File: node.h
// Author: Stanley Goodwin
// Creation Date: 6/5/2022
// Last Modified: 7/4/2022
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
	Node& normalize();  // Normalizes the current vector
	Node proj(Node v);  // Projects current vector on vector v
	double det();       // Finds the magnetude of the vector
	double det2();      // Finds the square magnetude of the vector 
};


// Node functions
double dot_product(Node v1, Node v2);
Node cross_product(Node v1, Node v2);


#endif NODE_H