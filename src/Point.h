#ifndef POINT_H
#define POINT_H
#pragma once


#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "NodeMesh.h"

class Point
{
public:
	Point();

	Point(const int& index,
		const std::string& name,
		const std::vector<double>& coordinates,
		const double& lcar,
		const bool& discretization = true);

	~Point();

	int getIndex();

	std::string getName();

	double getX();

	double getY();

	std::vector<double> getCoordinates();

	double getLcar();

	bool getDiscretization();

	NodeMesh getPointNode();

	std::string getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void setX(const double& x);

	void setY(const double& y);

	void setLcar(const double& lcar);

	void setDiscretization(const bool& discretization);

	void addNodeToPoint(NodeMesh* node);

private:
	int index_;           					// Point index
	std::string name_;    					// Gmsh Physical entity name
	std::vector<double> coordinates_;   	// Coordinates vector (x,y)
	double lcar_; 							// Characteristic length Gmsh parameter
	bool discretization_; 					// Choose to discretize a point with a mesh node
	NodeMesh* pointNode_;     					// Defines the mesh node discretizing the point
};


///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

Point::Point() {}

Point::Point(const int& index,
				const std::string& name,
				const std::vector<double>& coordinates,
				const double& lcar,
				const bool& discretization){
	index_ = index;
	name_ = name;
	coordinates_.reserve(2);
	coordinates_ = coordinates;
	lcar_ = lcar;
	discretization_ = discretization;
}

Point::~Point() {}

int Point::getIndex(){
	return index_;
}

std::string Point::getName(){
	return name_;
}

double Point::getX(){
	return coordinates_[0];
}

double Point::getY(){
	return coordinates_[1];
}

double Point::getLcar(){
	return lcar_;
}

bool Point::getDiscretization(){
	return discretization_;
}

NodeMesh Point::getPointNode(){
	return *pointNode_;
}

std::string Point::getGmshCode(){
	std::stringstream text;
	text << std::fixed;
	if (discretization_) {
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_
			<< "}; Physical Point('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", 0.0, " << lcar_ << "};\n//\n";
		return text.str();
	}
}

void Point::setIndex(const int& index){
	index_ = index;
}

void Point::setName(const std::string& name){
	name_ = name;
}

void Point::setX(const double& x){
	coordinates_[0] = x;
}

void Point::setY(const double& y){
	coordinates_[1] = y;
}

void Point::setLcar(const double& lcar){
	lcar_ = lcar;
}

void Point::setDiscretization(const bool& discretization){
	discretization_ = discretization;
}

void Point::addNodeToPoint(NodeMesh* node){
	pointNode_ = node;
}


#endif
