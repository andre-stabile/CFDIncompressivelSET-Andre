#pragma once

#include "LineLoop.h"
#include "Material.h"

class PlaneSurface
{
public:
	PlaneSurface();

	PlaneSurface(const int& index, const std::string& name, std::vector<LineLoop*> lineLoop, const double& thickness = 1.0);

	~PlaneSurface();

	int getIndex();

	std::string getName();

	std::vector<LineLoop*> getLineLoop();

	Material* getMaterial();

	double getThickness();

	std::vector<ElementMesh*> getElements();

	void setMaterial(Material* material);

	void setThickness(const double& thickness);

	void addElementToSurface(ElementMesh* object);

	std::string getGmshCode();

private:
	int index_;
	std::string name_;
	double thickness_;
	std::vector<LineLoop*> lineLoop_;
	Material* material_;
	std::vector<ElementMesh*> elements_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

PlaneSurface::PlaneSurface() {}

PlaneSurface::PlaneSurface(const int& index, const std::string& name, std::vector<LineLoop*> lineLoop, const double& thickness)
{
	index_ = index;
	name_ = name;
	lineLoop_ = lineLoop;
	thickness_ = thickness;
}

PlaneSurface::~PlaneSurface() {}

int PlaneSurface::getIndex()
{
	return index_;
}

std::string PlaneSurface::getName()
{
	return name_;
}

std::vector<LineLoop*> PlaneSurface::getLineLoop()
{
	return lineLoop_;
}

Material* PlaneSurface::getMaterial()
{
	return material_;
}

double PlaneSurface::getThickness()
{
	return thickness_;
}

std::vector<ElementMesh*> PlaneSurface::getElements()
{
	return elements_;
}

void PlaneSurface::setMaterial(Material* material)
{
	material_ = material;
}

void PlaneSurface::setThickness(const double& thickness)
{
	thickness_ = thickness;
}

void PlaneSurface::addElementToSurface(ElementMesh* object)
{
	elements_.push_back(object);
}

std::string PlaneSurface::getGmshCode()
{
	std::stringstream text;
	if (lineLoop_.size() == 1){
		text << name_ << " = news; Plane Surface(" << name_ << ") = {" << lineLoop_[0]->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	}else{
		text << name_ << " = news; Plane Surface(" << name_ << ") = {" ;
		for (int i=0; i<lineLoop_.size(); i++){
		 text << lineLoop_[i]->getName();
		 if(i!=(lineLoop_.size()-1)) text << " , " ;
		};
		text << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";

	};
	return text.str();
}