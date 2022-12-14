#ifndef NODE3DOF_H
#define NODE3DOF_H

# include <vector>

using namespace std;

class Node3DOF{
    public:
        Node3DOF(int index, const std::vector<double> &coords, const std::vector<double> &vel = {0.0, 0.0}, const std::vector<double> &accel = {0.0, 0.0});
        Node3DOF() {};
        ~Node3DOF();

        int getIndex();
        void setIndex(int index);

        std::vector<double> getInitialPositions();
        double getInitialPosition(int dir);
        void setInitialPositions(const std::vector<double> &coords);
        void setInitialPosition(int dir, double coord);

        std::vector<double> getCurrentPositions();
        double getCurrentPosition(int dir);
        void setCurrentPositions(const std::vector<double> &coords);
        void setCurrentPosition(int dir, double coord);

        std::vector<double> getVelocities();
        double getVelocity(int dir);
        void setVelocities(const std::vector<double> &vel);
        void setVelocity(int dir, double vel);

        std::vector<double> getAccelerations();
        double getAcceleration(int dir);
        void setAccelerations(const std::vector<double> &accel);
        void setAcceleration(int dir, double accel);

        void operator = (const Node3DOF &N);

    private:
        int index_;
        std::vector<double> initialPos_, currentPos_, vel_, accel_;
};

Node3DOF::Node3DOF(int index, const std::vector<double> &coords, const std::vector<double> &vel, const std::vector<double> &accel){
    index_ = index;
    initialPos_ = coords;
    currentPos_ = coords;
    vel_ = vel;
    accel_ = accel;
}

Node3DOF::~Node3DOF() {}

int Node3DOF::getIndex(){
    return index_;
}
void  Node3DOF::setIndex(int index){
    index_ = index;
}

std::vector<double> Node3DOF::getInitialPositions(){
    return initialPos_;
}
double Node3DOF::getInitialPosition(int dir){
    return initialPos_[dir];
}
void Node3DOF::setInitialPositions(const std::vector<double> &coords){
    initialPos_ = coords;
}
void Node3DOF::setInitialPosition(int dir, double coord){
    initialPos_[dir] = coord;
}

std::vector<double> Node3DOF::getCurrentPositions(){
    return currentPos_;
}
double Node3DOF::getCurrentPosition(int dir){
    return currentPos_[dir];
}
void Node3DOF::setCurrentPositions(const std::vector<double> &coords){
    currentPos_ = coords;
}
void Node3DOF::setCurrentPosition(int dir, double coord){
    currentPos_[dir] = coord;
}

std::vector<double> Node3DOF::getVelocities(){
    return vel_;
}
double Node3DOF::getVelocity(int dir){
    return vel_[dir];
}
void Node3DOF::setVelocities(const std::vector<double> &vel){
    vel_ = vel;
}
void Node3DOF::setVelocity(int dir, double vel){
    vel_[dir] = vel;
}

std::vector<double> Node3DOF::getAccelerations(){
    return accel_;
}
double Node3DOF::getAcceleration(int dir){
    return accel_[dir];
}
void Node3DOF::setAccelerations(const std::vector<double> &accel){
    accel_ = accel;
}
void Node3DOF::setAcceleration(int dir, double accel){
    accel_[dir] = accel;
}

void Node3DOF::operator = (const Node3DOF &N){
    index_ = N.index_;
    initialPos_ = N.initialPos_;
    currentPos_ = N.initialPos_;
    vel_ = N.vel_;
    accel_ = N.accel_;
}

#endif