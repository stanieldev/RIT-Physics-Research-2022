/*
 * File:	node.hpp
 * Author:	Stanley Goodwin
 * 
 * Stores all the necessary data of the 3-component vector, 
 * Allows operations between node-node and node-scalar.
 * Includes some other vector functions like projection, normalization, and magnetude.
 */
#pragma once
#ifndef NODE_H
#define NODE_H

struct Node {
	double x, y, z;

	Node();
	Node(const Node& _node);
	Node(double _x, double _y, double _z);
	
	Node& operator= (const Node& _node);
	Node  operator+ (const Node& _node);
	Node& operator+=(const Node& _node);
	Node  operator- (const Node& _node);
	Node& operator-=(const Node& _node);
	Node  operator* (double _scalar);
	Node& operator*=(double _scalar);
	Node  operator/ (double _scalar);
	Node& operator/=(double _scalar);
	
	Node normalize();
	Node project(Node _onto);
	double magnetude();
	double magnetude_squared();
};

Node normalize(Node _node);
Node project(Node _node, Node _onto);
double dot_product(Node _node1, Node _node2);
Node cross_product(Node _node1, Node _node2);

#endif NODE_H