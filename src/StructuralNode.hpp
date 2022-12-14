#ifndef STRUCTURALNODE_H
#define STRUCTURALNODE_H

# include <vector>

using namespace std;

class StructuralNode{
    public:
        StructuralNode(int index, const std::vector<double> &coords, const std::vector<double> &vel = {0.0, 0.0}, const std::vector<double> &accel = {0.0, 0.0});

        ~StructuralNode();

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

        void operator = (const StructuralNode &N);

    private:
        int index_;
        std::vector<double> initialPos_, currentPos_, vel_, accel_;
};

StructuralNode::StructuralNode(int index, const std::vector<double> &coords, const std::vector<double> &vel, const std::vector<double> &accel){
    index_ = index;
    initialPos_ = coords;
    currentPos_ = coords;
    vel_ = vel;
    accel_ = accel;
}

StructuralNode::~StructuralNode() {}

int StructuralNode::getIndex(){
    return index_;
}
void  StructuralNode::setIndex(int index){
    index_ = index;
}

std::vector<double> StructuralNode::getInitialPositions(){
    return initialPos_;
}
double StructuralNode::getInitialPosition(int dir){
    return initialPos_[dir];
}
void StructuralNode::setInitialPositions(const std::vector<double> &coords){
    initialPos_ = coords;
}
void StructuralNode::setInitialPosition(int dir, double coord){
    initialPos_[dir] = coord;
}

std::vector<double> StructuralNode::getCurrentPositions(){
    return currentPos_;
}
double StructuralNode::getCurrentPosition(int dir){
    return currentPos_[dir];
}
void StructuralNode::setCurrentPositions(const std::vector<double> &coords){
    currentPos_ = coords;
}
void StructuralNode::setCurrentPosition(int dir, double coord){
    currentPos_[dir] = coord;
}

std::vector<double> StructuralNode::getVelocities(){
    return vel_;
}
double StructuralNode::getVelocity(int dir){
    return vel_[dir];
}
void StructuralNode::setVelocities(const std::vector<double> &vel){
    vel_ = vel;
}
void StructuralNode::setVelocity(int dir, double vel){
    vel_[dir] = vel;
}

std::vector<double> StructuralNode::getAccelerations(){
    return accel_;
}
double StructuralNode::getAcceleration(int dir){
    return accel_[dir];
}
void StructuralNode::setAccelerations(const std::vector<double> &accel){
    accel_ = accel;
}
void StructuralNode::setAcceleration(int dir, double accel){
    accel_[dir] = accel;
}

void StructuralNode::operator = (const StructuralNode &N){
    index_ = N.index_;
    initialPos_ = N.initialPos_;
    currentPos_ = N.initialPos_;
    vel_ = N.vel_;
    accel_ = N.accel_;
}

#endif