/*
 * File:	node.hpp
 * Author:	Stanley Goodwin
 * 
 * Stores all the necessary data of the 3-component
 * vector, as well as the operations between node-node
 * and node-scalar. Includes some other vector functions
 * like projection, normalization, and determinant.
 */
#pragma once
#ifndef NODE_H
#define NODE_H

struct Node {
	double x, y, z;

	Node();
	Node(double ix, double iy, double iz);
	Node(const Node& node);

	Node& operator= (const Node& node);
	Node  operator+ (const Node& node);
	Node& operator+=(const Node& node);
	Node  operator- (const Node& node);
	Node& operator-=(const Node& node);
	Node  operator* (double scalar);
	Node& operator*=(double scalar);
	Node  operator/ (double scalar);
	Node& operator/=(double scalar);

	Node& normalize();           // Normalizes the current node
	void project(Node onto);     // Projects current vector on node v
	double magnetude();          // Finds the magnetude of the node
	double magnetude_squared();  // Finds the square magnetude of the node
};

double dot_product(Node node1, Node node2);
Node cross_product(Node node1, Node node2);
Node project(Node node, Node onto);

#endif NODE_H