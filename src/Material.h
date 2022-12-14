#pragma once

#include <string>

class Material
{
public:
	Material(const int& index, const double& young, const double& poisson, const double& density = 0.0, const std::string& behaviour = "PLANE_STRESS");

	~Material();

	int getIndex();

	double getYoung();

	double getPoisson();

	double getDensity();

	std::string getBehaviour();

	void setIndex(const int& index);
	
	void setYoung(const double& young);

	void setPoisson(const double& poisson);

	void setDensity(const double& density);

	void setBehaviour(const std::string& behaviour);

private:
	int index_;
	double young_;
	double poisson_;
	double density_;
	std::string behaviour_;
};

///----------------------------------------------------------------------------
///-------------------------------IMPLEMENTATION-------------------------------
///----------------------------------------------------------------------------

Material::Material(const int& index, const double& young, const double& poisson, const double& density, const std::string& behaviour)
{
	index_ = index;
	young_ = young;
	poisson_ = poisson;
	density_ = density;
	behaviour_ = behaviour;
}

Material::~Material() {}

int Material::getIndex()
{
	return index_;
}

double Material::getYoung()
{
	return young_;
}

double Material::getPoisson()
{
	return poisson_;
}

double Material::getDensity()
{
	return density_;
}

std::string Material::getBehaviour()
{
	return behaviour_;
}

void Material::setIndex(const int& index)
{
	index_ = index;
}

void Material::setYoung(const double& young)
{
	young_ = young;
}

void Material::setPoisson(const double& poisson)
{
	poisson_ = poisson;
}

void Material::setDensity(const double& density)
{
	density_ = density;
}

void Material::setBehaviour(const std::string& behaviour)
{
	behaviour_ = behaviour;
}