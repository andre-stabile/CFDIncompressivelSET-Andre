#pragma once

#include "Point.h"
#include "ElementMesh.h"

class Line
{
public:
	Line();

	Line(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization = true);

	~Line();

	Line* operator-();

	int getIndex();

	std::string getName();

	Point* getInitialPoint();

	Point* getEndPoint();

	bool getDiscretization();

	std::vector<NodeMesh*> getLineNodes();

	std::string getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void setInitialPoint(Point& point);

	void setEndPoint(Point& point);

	void setDiscretization(const bool& discretization);

	void addNodesToLine(const std::vector<NodeMesh*>& nodes);

	void addElementsToLine(ElementMesh* element);

private:
	int index_;
	std::string name_;
	std::vector<Point*> points_;
	bool discretization_;
	std::vector<NodeMesh*> lineNodes_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

Line::Line() {}

Line::Line(const int& index, const std::string& name, std::vector<Point*> points, const bool& discretization)
{
	index_ = index;
	name_ = name;
	discretization_ = discretization;
	points_.reserve(2);
		for (Point* point : points)
			points_.push_back(point);
}

Line::~Line() {}

Line* Line::operator-()
{
//    if(name_[0] == '-') {
//	Line* copy = this;
//	copy->setName(name_.erase(0));
//	return *copy;
//    } else {
//	Line* copy = this;
//	copy->setName(name_.insert(0, "-"));
//	return *copy;
//    }
	Line* copy = new Line(index_, name_, {points_});
	copy->setName("-" + name_);
	return copy;
}

int Line::getIndex()
{
	return index_;
}
std::string Line::getName()
{
	return name_;
}
Point* Line::getInitialPoint()
{
	return points_[0];
}
Point* Line::getEndPoint()
{
	return points_[1];
}
bool Line::getDiscretization()
{
	return discretization_;
}
std::vector<NodeMesh*> Line::getLineNodes()
{
	return lineNodes_;
}
std::string Line::getGmshCode()
{
	std::stringstream text;

	if (discretization_) {
		if(points_.size() == 3){
			text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName() << ", " << points_[2]->getName()
				<< "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		}else{
			text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() ;
			for (int i = 1; i < points_.size(); ++i)
			{
				text <<  ", " << points_[i]->getName();	
			}
			text << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		}
		return text.str();
	}
	else {
		if(points_.size() == 3){
			text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName() << ", " << points_[2]->getName()
				<< "};\n//\n";
		}else{
		text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() ;
			for (int i = 1; i < points_.size(); ++i)
			{
				text <<  ", " << points_[i]->getName();	
			}
			text << "}; ";
		}
		return text.str();
	}
}
void Line::setIndex(const int& index)
{
	index_ = index;
}
void Line::setName(const std::string& name)
{
	name_ = name;
}
void Line::setInitialPoint(Point& point)
{
	//points_[0] = point;
}
void Line::setEndPoint(Point& point)
{
	//points_[1] = point;
}
void Line::setDiscretization(const bool& discretization)
{
	discretization_ = discretization;
}
void Line::addNodesToLine(const std::vector<NodeMesh*>& nodes)
{
	for (NodeMesh* node : nodes)
	lineNodes_.push_back(node);
}
void Line::addElementsToLine(ElementMesh* element)
{

}