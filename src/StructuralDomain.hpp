#ifndef STRUCTURALDOMAIN_H
#define STRUCTURALDOMAIN_H

#include "StructuralNode.hpp"
#include "Node3DOF.hpp"
#include <cmath>

using namespace std;

class StructuralDomain{
    public:
        StructuralDomain(const Node3DOF &center, const std::vector<StructuralNode> &boundary, const std::vector<bool> &constraints,
                        const std::vector<double> &springs, const std::vector<double> &mass, const std::vector<double> &dampers,
                        double dt = 0.0, double beta = 0.25, double gamma = 0.5);
        StructuralDomain() {};
        ~StructuralDomain();

        Node3DOF getCenter();
        void setCenter(const Node3DOF &center);

        std::vector<StructuralNode> getBoundary();
        void setBoundary(const std::vector<StructuralNode> &boundary);

        std::vector<double> getMasses();
        double getMass(int dir);
        void setMasses(const std::vector<double> &mass);
        void setMass(int dir, double mass);

        std::vector<double> getDampers();
        double getDamper(int dir);
        void setDampers(const std::vector<double> &dampers);
        void setDamper(int dir, double dampers);

        std::vector<bool> getConstraints();
        bool getConstraint(int dir);
        void setConstraints(const std::vector<bool> &constraints);
        void setConstraint(int dir, bool constraint);

        std::vector<double> getSprings();
        double getSpring(int dir);
        void setSprings(const std::vector<double> &springs);
        void setSpring(int dir, double spring);

        double getDt();
        void setDt(double dt);

        double getBeta();
        void setBeta(double beta);

        double getGamma();
        void setGamma(double gamma);

        std::vector<double> getRadii();
        std::vector<double> getAngles();

        void updateCenterPosition(const std::vector<double> &f);
        void updateBoundary();

    private:
        double dt_, beta_, gamma_;
        Node3DOF center_;
        std::vector<StructuralNode> boundaryNodes_;
        std::vector<bool> constraints_;
        std::vector<double> springs_, mass_, dampers_, radii_, angles_;
};

StructuralDomain::StructuralDomain(const Node3DOF &center, const std::vector<StructuralNode> &boundary, const std::vector<bool> &constraints,
                                   const std::vector<double> &springs, const std::vector<double> &mass, const std::vector<double> &dampers,
                                   double dt, double beta, double gamma)
{
    center_ = center;
    boundaryNodes_ = boundary;
    constraints_ = constraints;

    springs_ = springs;
    mass_ = mass;
    dampers_ = dampers;


    for (int i = 0; i < boundaryNodes_.size(); i++){
        radii_.push_back(sqrt(pow(boundaryNodes_[i].getInitialPosition(0)-center_.getInitialPosition(0),2)+pow(boundaryNodes_[i].getInitialPosition(1)-center_.getInitialPosition(1),2)));
        angles_.push_back(atan2(boundaryNodes_[i].getInitialPosition(1)-center_.getInitialPosition(1),boundaryNodes_[i].getInitialPosition(0)-center_.getInitialPosition(0)));
    }

    dt_ = dt;
    beta_ = beta;
    gamma_ = gamma;
}

StructuralDomain::~StructuralDomain() {}

Node3DOF StructuralDomain::getCenter() {
    return center_;
}
void StructuralDomain::setCenter(const Node3DOF &center){
    center_ = center;
}

std::vector<StructuralNode> StructuralDomain::getBoundary(){
    return boundaryNodes_;
}
void StructuralDomain::setBoundary(const std::vector<StructuralNode> &boundary){
    boundaryNodes_ = boundary;
}

std::vector<double> StructuralDomain::getMasses(){
    return mass_;
}
double StructuralDomain::getMass(int dir){
    return mass_[dir];
}
void StructuralDomain::setMasses(const std::vector<double> &mass){
    mass_ = mass;
}
void StructuralDomain::setMass(int dir, double mass){
    mass_[dir] = mass;
}

std::vector<double> StructuralDomain::getDampers(){
    return dampers_;
}
double StructuralDomain::getDamper(int dir){
    return dampers_[dir];
}
void StructuralDomain::setDampers(const std::vector<double> &dampers){
    dampers_ = dampers;
}
void StructuralDomain::setDamper(int dir, double damper){
    dampers_[dir] = damper;
}

std::vector<bool> StructuralDomain::getConstraints(){
    return constraints_;
}
bool StructuralDomain::getConstraint(int dir){
    return constraints_[dir];
}
void StructuralDomain::setConstraints(const std::vector<bool> &constraints){
    constraints_ = constraints;
}
void StructuralDomain::setConstraint(int dir, bool constraint){
    constraints_[dir] = constraint;
}

std::vector<double> StructuralDomain::getSprings(){
    return springs_;
}
double StructuralDomain::getSpring(int dir){
    return springs_[dir];
}
void StructuralDomain::setSprings(const std::vector<double> &springs){
    springs_ = springs;
}
void StructuralDomain::setSpring(int dir, double spring){
    springs_[dir] = spring;
}

double StructuralDomain::getDt(){
    return dt_;
}
void StructuralDomain::setDt(double dt){
    dt_ = dt;
}

double StructuralDomain::getBeta(){
    return beta_;
}
void StructuralDomain::setBeta(double beta){
   beta_ = beta;
}

double StructuralDomain::getGamma(){
    return gamma_;
}
void StructuralDomain::setGamma(double gamma){
    gamma_ = gamma;
}

std::vector<double> StructuralDomain::getRadii(){
    return radii_;
}
std::vector<double> StructuralDomain::getAngles(){
    return angles_;
}

void StructuralDomain::updateCenterPosition(const std::vector<double> &f){
    for (int i = 0; i < 3; i++){
        if (!constraints_[i]){
            double y = center_.getCurrentPosition(i);
            double x = center_.getInitialPosition(i);
            double v = center_.getVelocity(i);
            double a = center_.getAcceleration(i);
            double k = springs_[i];
            double c = dampers_[i];
            double m = mass_[i];

            std::vector< std::vector <double> > A(2);
            A[0].resize(2);
            A[1].resize(2);
            A[0][0] = dt_*dt_*beta_*k/m + 1.0;
            A[0][1] = dt_*dt_*beta_*c/m;
            A[1][0] = dt_*gamma_*k/m;
            A[1][1] = dt_*gamma_*c/m + 1.0;
            double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

            std::vector<double> b(2);
            std::vector<double> res(2);
            b[0] = y + dt_*v + dt_*dt_*((0.5-beta_)*a + beta_*(f[i]+k*x)/m);
            b[1] = v + dt_*((1-gamma_)*a + gamma_*(f[i]+k*x)/m);

            res[0] = (A[1][1]*b[0] - A[0][1]*b[1])/detA;
            res[1] = (A[0][0]*b[1] - A[1][0]*b[0])/detA;

            a = f[i]/m - c*res[1]/m - k*(res[0]-x)/m;

            center_.setCurrentPosition(i, res[0]);
            center_.setVelocity(i, res[1]);
            center_.setAcceleration(i, a);
        }
    }       
}

void StructuralDomain::updateBoundary(){
    for (int i = 0; i < boundaryNodes_.size(); i++){
        boundaryNodes_[i].setCurrentPosition(0, center_.getCurrentPosition(0) + radii_[i] * cos(angles_[i]+center_.getCurrentPosition(2)));
        boundaryNodes_[i].setCurrentPosition(1, center_.getCurrentPosition(1) + radii_[i] * sin(angles_[i]+center_.getCurrentPosition(2)));
        boundaryNodes_[i].setVelocity(0, center_.getVelocity(0) - radii_[i] * sin(angles_[i]+center_.getCurrentPosition(2)) * center_.getVelocity(2));
        boundaryNodes_[i].setVelocity(1, center_.getVelocity(1) + radii_[i] * cos(angles_[i]+center_.getCurrentPosition(2)) * center_.getVelocity(2));
        boundaryNodes_[i].setAcceleration(0, center_.getAcceleration(0) - radii_[i] * (cos(angles_[i]+center_.getCurrentPosition(2)) * pow(center_.getVelocity(2),2) +
                                                                                       sin(angles_[i]+center_.getCurrentPosition(2)) * center_.getAcceleration(2)));
        boundaryNodes_[i].setAcceleration(1, center_.getAcceleration(1) + radii_[i] * (cos(angles_[i]+center_.getCurrentPosition(2)) * center_.getAcceleration(2) -
                                                                                       sin(angles_[i]+center_.getCurrentPosition(2)) * pow(center_.getVelocity(2),2)));
    }
}

#endif