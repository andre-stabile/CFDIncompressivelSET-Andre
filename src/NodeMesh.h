#pragma once

#include <vector>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric::ublas;

class NodeMesh
{
public:
	NodeMesh(const int& index, const bounded_vector<double,2>& coord);

	~NodeMesh();

	int getIndex();

	bounded_vector<double, 2> getCoordinates();

	double getX();

	double getY();

	void setIndex(const int& index);

	void setCoordinates(const bounded_vector<double, 2>& coord);

	void setX(const double& x);

	void setY(const double& y);

private:
	int index_;
	bounded_vector<double, 2> coord_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

NodeMesh::NodeMesh(const int& index, const bounded_vector<double,2>& coord)
{
	index_ = index;
	coord_ = coord;
}

NodeMesh::~NodeMesh() {}

int NodeMesh::getIndex()
{
	return index_;
}

double NodeMesh::getX()
{
	return coord_(0);
}

double NodeMesh::getY()
{
	return coord_(1);
}

bounded_vector<double, 2> NodeMesh::getCoordinates()
{
	return coord_;
}

void NodeMesh::setIndex(const int& index)
{
	index_ = index;
}

void NodeMesh::setX(const double& x)
{
	coord_(0) = x;
}

void NodeMesh::setY(const double& y)
{
	coord_(1) = y;
}