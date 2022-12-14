#pragma once

#include "PlaneSurface.h"
#include "BoundaryCondition.h"
#include <string>

class Geometry
{
public:
	Geometry();

	Geometry(const int& index);

	~Geometry();

	int getIndex();

	int getNumberOfPoints();

	int getNumberOfLines();

	int getNumberOfLineLoops();

	int getNumberOfPlaneSurfaces();

	int getNumberOfBoundaryConditions(const std::string& type);

	Point* getPoint(const std::string& name);

	Line* getLine(const std::string& name);

	LineLoop* getLineLoop(const std::string& name);

	PlaneSurface* getPlaneSurface(const std::string& name);

	std::unordered_map<std::string, PlaneSurface*> getPlaneSurfaces();

	std::string getGmshCode();

	Point* addPoint(std::vector<double> coordinates, const double& lcar = 1.0, const bool& discretization = true);

	Line* addLine(std::vector<Point*> points, const bool& discretization = true);

	LineLoop* addLineLoop(std::vector<Line*> lines);

	PlaneSurface* addPlaneSurface(std::vector<LineLoop*> lineLoop, double thickness = 1.0);

	PlaneSurface* addPlaneSurface(std::vector<Line*> lines, double thickness = 1.0);

	void appendGmshCode(std::string text);

	void transfiniteLine(std::vector<Line*> lines, const int& divisions, const double& progression = 1);

	void transfiniteSurface(std::vector<PlaneSurface*> surfaces, std::string oientation = "Left", std::vector<Point*> points = std::vector<Point*>());

	void addBoundaryCondition(const std::string& type, Point* point, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);
	
	void addBoundaryCondition(const std::string& type, Line* line, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);

	std::vector<BoundaryCondition*> getBoundaryCondition(const std::string& type);

private:
	int index_;
	std::unordered_map<std::string, Point*> points_;
	std::unordered_map<std::string, Line*> lines_;
	std::unordered_map<std::string, LineLoop*> lineLoops_;
	std::unordered_map<std::string, PlaneSurface*> planeSurfaces_;
	std::unordered_map<std::string, std::vector<BoundaryCondition*>> boundaryConditions_;
	std::string gmshCode_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

Geometry::Geometry() {}

Geometry::Geometry(const int& index)
{
	index_ = index;
}

Geometry::~Geometry() {}

int Geometry::getIndex()
{
	return index_;
}

int Geometry::getNumberOfPoints()
{
	return points_.size();
}

int Geometry::getNumberOfLines()
{
	return lines_.size();
}

int Geometry::getNumberOfLineLoops()
{
	return lineLoops_.size();
}

int Geometry::getNumberOfPlaneSurfaces()
{
	return planeSurfaces_.size();
}

int Geometry::getNumberOfBoundaryConditions(const std::string& type)
{
	return boundaryConditions_[type].size();
}


Point* Geometry::getPoint(const std::string& name)
{
	return points_[name];
}

Line* Geometry::getLine(const std::string& name)
{
	return lines_[name];
}

LineLoop* Geometry::getLineLoop(const std::string& name)
{
	return lineLoops_[name];
}

PlaneSurface* Geometry::getPlaneSurface(const std::string& name)
{
	return planeSurfaces_[name];
}

std::unordered_map<std::string, PlaneSurface*> Geometry::getPlaneSurfaces()
{
	return planeSurfaces_;
}

std::string Geometry::getGmshCode()
{
	return gmshCode_;
}

Point* Geometry::addPoint(std::vector<double> coordinates, const double& lcar, const bool& discretization)
{
	int index = getNumberOfPoints();
	std::stringstream name;
	name << "p" << index;
	Point* p = new Point(index, name.str(), coordinates, lcar, discretization);
	points_[p->getName()] = p;
	gmshCode_ += p->getGmshCode();
	return p;
}

Line* Geometry::addLine(std::vector<Point*> points, const bool& discretization)
{
	int index = getNumberOfLines();
	std::stringstream name;
	name << "l" << index;
	Line* l = new Line(index, name.str(), points, discretization);
	lines_[l->getName()] = l;
	gmshCode_ += l->getGmshCode();
	return l;
}

LineLoop* Geometry::addLineLoop(std::vector<Line*> lines)
{
	int index = getNumberOfLineLoops();
	std::stringstream name;
	name << "ll" << index;
	LineLoop* ll = new LineLoop(index, name.str(), lines);
	ll->verification();
	lineLoops_[ll->getName()] = ll;
	gmshCode_ += ll->getGmshCode();
	return ll;
}

PlaneSurface* Geometry::addPlaneSurface(std::vector<LineLoop*> lineLoop, double thickness)
{
	int index = getNumberOfPlaneSurfaces();
	std::stringstream name;
	name << "s" << index;
	PlaneSurface* s = new PlaneSurface(index, name.str(), lineLoop, thickness);
	planeSurfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

PlaneSurface* Geometry::addPlaneSurface(std::vector<Line*> lines, double thickness)
{
	int index = getNumberOfPlaneSurfaces();
	std::stringstream name;
	name << "s" << index;
	LineLoop* ll = addLineLoop(lines);
	PlaneSurface* s = new PlaneSurface(index, name.str(), {ll}, thickness);
	planeSurfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

void Geometry::appendGmshCode(std::string text)
{
	gmshCode_ += text;
}

void Geometry::transfiniteLine(std::vector<Line*> lines, const int& divisions, const double& progression)
{
	std::stringstream text;
	text << "Transfinite Line {";
	for (size_t i = 0; i < lines.size(); i++)
	{
		text << lines[i]->getName();
		if (i != (lines.size() - 1))
			text << ", ";
	}
	text << "} = " << divisions << " Using Progression " << progression << ";\n//\n";
	gmshCode_ += text.str();
}

void Geometry::transfiniteSurface(std::vector<PlaneSurface*> planeSurfaces, std::string orientation, std::vector<Point*> points)
{
	std::stringstream text;
	text << "Transfinite Surface {";
	for (size_t i = 0; i < planeSurfaces.size(); i++)
	{
		text << planeSurfaces[i]->getName();
		if (i != (planeSurfaces.size() - 1))
			text << ", ";
	}
	text << "} ";
	if (points.size() != 0)
	{
		text << "= {";
		for (size_t i = 0; i < points.size(); i++)
		{
			text << points[i]->getName();
			if (i != (points.size() - 1))
				text << ", ";
		}
		text << "} " << orientation << ";\n//\n";
	}
	else
	{
		text << orientation << ";\n//\n";
	}
	gmshCode_ += text.str();
}

void Geometry::addBoundaryCondition(const std::string& type, Point* point, const std::vector<double>& componentX, const std::vector<double>& componentY, const std::string& referenceSystem, const std::string& method, const double& penaltyParameter)
{
	if (boundaryConditions_.count(type) == 0)
		boundaryConditions_[type] = std::vector<BoundaryCondition*>();
	
	int index = getNumberOfBoundaryConditions(type);
	BoundaryCondition* b = new BoundaryCondition(index, point->getName(), componentX, componentY, referenceSystem, method, penaltyParameter);
	boundaryConditions_[type].push_back(b);
}

void Geometry::addBoundaryCondition(const std::string& type, Line* line, const std::vector<double>& componentX, const std::vector<double>& componentY, const std::string& referenceSystem, const std::string& method, const double& penaltyParameter)
{
	if (boundaryConditions_.count(type) == 0)
		boundaryConditions_[type] = std::vector<BoundaryCondition*>();

	int index = getNumberOfBoundaryConditions(type);
	BoundaryCondition* b = new BoundaryCondition(index, line->getName(), componentX, componentY, referenceSystem, method, penaltyParameter);
	boundaryConditions_[type].push_back(b);
}


std::vector<BoundaryCondition*> Geometry::getBoundaryCondition(const std::string& type)
{
	return boundaryConditions_[type];
}