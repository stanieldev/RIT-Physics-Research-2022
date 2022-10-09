/**
 * @file	node.h
 * @brief	A class of 3-vectors and operations thereof.
 *
 * @author	Stanley Goodwin
 * Contact: sfg99709akwork@gmail.com
 *
 * Creation Date: 6/5/2022
 * Last Modified: 7/18/2022
 */
#pragma once
//#ifndef NODE_H
#define NODE_H


/**
 * Node characteristics class.
 *
 * Stores all the necessary data of the 3-component
 * vector, as well as the operations between node-node
 * and node-scalar. Includes some other vector functions
 * like projection, normalization, and determinant.
 *
 */
class Node {
public:
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
	Node& normalize();  // Normalizes the current node
	Node proj(Node node);  // Projects current vector on node v
	double det();       // Finds the magnetude of the node
	double det2();      // Finds the square magnetude of the node
};


// Extra vector operations
double dot_product(Node node1, Node node2);  // Vector dot product
Node cross_product(Node node1, Node node2);  // Vector cross product


//#endif NODE_H