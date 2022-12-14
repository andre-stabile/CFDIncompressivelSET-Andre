#pragma once

#include <string>
#include <vector>

class BoundaryCondition
{
public:
	BoundaryCondition(const int& index, const std::string& object, const std::vector<double>& componentX = std::vector<double>(), const std::vector<double>& componentY = std::vector<double>(), const std::string& referenceSystem = "GLOBAL", const std::string& method = "STRONG", const double& penaltyParameter = 1.0e6);

	~BoundaryCondition();

	int getIndex();

	std::string getPointName();

	std::string getLineName();

	std::string getReferenceSystem();

	std::vector<double> getComponentX();

	std::vector<double> getComponentY();

	std::string getMethod();

	double getPenaltyParameter();

	void setIndex(const int& index);

	void setPointName(const std::string& name);

	void setLineName(const std::string& name);

	void setReferenceSystem(const std::string& referenceSystem);

	void setComponentX(std::vector<double> componentX);

	void setComponentY(std::vector<double> componentY);

	void setMethod(const std::string& methhod);

	void setPenaltyParameter(const double& penaltyParameter);

private:
	int index_;
	std::string point_;
	std::string line_;
	std::string referenceSystem_;
	std::vector<double> componentX_;
	std::vector<double> componentY_;
	std::string method_;
	double penaltyParameter_;

};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

BoundaryCondition::BoundaryCondition(const int& index,
	const std::string& object,
	const std::vector<double>& componentX,
	const std::vector<double>& componentY,
	const std::string& referenceSystem,
	const std::string& method,
	const double& penaltyParameter)
{
	index_ = index;
	(object[0] == 'p') ? point_ = object : line_ = object;
	componentX_ = componentX;
	componentY_ = componentY;
	referenceSystem_ = referenceSystem;
	method_ = method;
	penaltyParameter_ = penaltyParameter;
}

BoundaryCondition::~BoundaryCondition() {}

int BoundaryCondition::getIndex()
{
	return index_;
}

std::string BoundaryCondition::getPointName()
{
	return point_;
}

std::string BoundaryCondition::getLineName()
{
	return line_;
}

std::string BoundaryCondition::getReferenceSystem()
{
	return referenceSystem_;
}

std::vector<double> BoundaryCondition::getComponentX()
{
	return componentX_;
}

std::vector<double> BoundaryCondition::getComponentY()
{
	return componentY_;
}

std::string BoundaryCondition::getMethod()
{
	return method_;
}

double BoundaryCondition::getPenaltyParameter()
{
	return penaltyParameter_;
}

void BoundaryCondition::setIndex(const int& index)
{
	index_ = index;
}

void BoundaryCondition::setPointName(const std::string& name)
{
	point_ = name;
}

void BoundaryCondition::setLineName(const std::string& name)
{
	line_ = name;
}

void BoundaryCondition::setReferenceSystem(const std::string& referenceSystem)
{
	referenceSystem_ = referenceSystem;
}

void BoundaryCondition::setComponentX(std::vector<double> componentX)
{
	componentX_ = componentX;
}

void BoundaryCondition::setComponentY(std::vector<double> componentY)
{
	componentY_ = componentY;
}

void BoundaryCondition::setMethod(const std::string& method)
{
	method_ = method;
}

void BoundaryCondition::setPenaltyParameter(const double& penaltyParameter)
{
	penaltyParameter_ = penaltyParameter;
}
