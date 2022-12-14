#pragma once

#include "NodeMesh.h"
#include "Material.h"

class ElementMesh
{
public:
	ElementMesh(const int& index, const std::vector<NodeMesh*>& nodes, Material* material, const double& thickness);

	~ElementMesh();

private:
	int index_;
	std::vector<NodeMesh*> nodes_;
	Material* material_;
	double thickness_;
	std::string elementType_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

ElementMesh::ElementMesh(const int& index, const std::vector<NodeMesh*>& nodes, Material* material, const double& thickness)
{
	index_ = index;
	nodes_ = nodes;
	material_ = material;
	thickness_ = thickness;
}

ElementMesh::~ElementMesh() {}
