//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------NODES-------------------------------------
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric;

/// Defines the node object and stores all nodal variables information

template<int DIM>

class Node{
public:
    /// Defines a type to store double variables with the dimension size
    typedef ublas::bounded_vector<double,DIM>        VecLocD;

    /// Defines a pair with integer and  double variables
    /// which helps to setting the constraints
    typedef std::vector<std::pair<int, double> >     VecConstrD;

    /// Defines the node iterator for helping to build loops in order
    /// to enforce bcs
    typedef VecConstrD::iterator                     VecConstrDIt;

    static const int spaceDim = DIM;

private:
    //Main variables
    VecLocD          coord_;                  //Nodal coordinate vector
    VecLocD          initialCoord_;           //Initial nodal coordinate vector
    VecLocD          coordUpdated_;           //Updated nodal coordinate vector
    int              index_;                  //Node index
    VecLocD          previousCoord_;          //Previous nodal coordinate vector
    VecLocD          nNodal_;                 //Nodal normal vector
    
    //Fluid
    int              constrainType[3];        //Constrain direction
    double           constrainValue[3];       //Nodal prescribed velocity

    int              constrainTypeLaplace[3]; //Constrain direction Laplace
    double           constrainValueLaplace[3];//Nodal prescribed value

    VecLocD          velocity_;               //Nodal velocity
    VecLocD          previousVelocity_;       //Previous time step velocity
    VecLocD          intermediateVelocity_;       //Intermediate velocity

    VecLocD          acceleration_;           //Nodal acceleration
    VecLocD          previousAcceleration_;   //Previous time step acceleration
    VecLocD          intermediateAcceleration_;//Intermediate acceleration

    double           pressure_;               //Nodal pressure 
    double           previousPressure_;       //Previous time step pressure

    double           vorticity_;              //Vorticity

    VecLocD          meshVelocity_;           //Nodal mesh velocity
    VecLocD          previousMeshVelocity_;   //Previous time step mesh velocity

    std::vector<int> invIncidence;
    
public:
    ///Constructor - Defines a node with index and coordinates
    Node(VecLocD& coor, int index){
        coord_ = coor;
        initialCoord_ = coor;
        coordUpdated_ = coor;
        previousCoord_ = coor;
        index_ = index; 

        constrainType[0] = 0;    constrainType[1] = 0;    constrainType[2] = 0;
        constrainValue[0] = 0;   constrainValue[1] = 0;   constrainValue[2] = 0;
        pressure_ = 0.;          previousPressure_ = 0.;  
        vorticity_ = 0.;
        velocity_.clear();     previousVelocity_.clear();
        intermediateVelocity_.clear();
        acceleration_.clear(); previousAcceleration_.clear(); 
        intermediateAcceleration_.clear();

        constrainTypeLaplace[0] = 0;    constrainTypeLaplace[1] = 0;
        constrainTypeLaplace[2] = 0;    constrainValueLaplace[0] = 0;
        constrainValueLaplace[1] = 0;   constrainValueLaplace[2] = 0;
        
        invIncidence.clear();

    };

    /// Clear all node object variables
    void clearVariables();

    /// Returns the node coordinate vector
    /// @return node coordinate vector
    VecLocD getCoordinates() {return coord_;};

    /// Returns the node coordinate component value
    /// @return node coordinate component value
    double getCoordinateValue(int dir) {return coord_(dir);};

    /// Returns the node initial coordinate vector
    /// @return node initial coordinate vector
    VecLocD getInitialCoordinates() {return initialCoord_;};

    /// Returns the node coordinate vector at the previous time step
    /// @return node coordinate vector at the previous time step
    VecLocD getPreviousCoordinates() {return previousCoord_;}

    /// Returns the node updated coordinate vector
    /// @return node coordinate updated vector
    VecLocD getUpdatedCoordinates() {return coordUpdated_;};

    /// Increment the coordinate vector
    /// @param int direction @param double increment value
    void incrementCoordinate(int dir, double u);

    /// Sets the previous coordinate vector
    /// @param int direction @param double value
    void setPreviousCoordinates(int dir, double u);

    /// Sets the node coordinate vector
    /// @param VecLocD Coordinate
    void setCoordinates(VecLocD& coor){coord_ = coor;};

    /// Sets the updated coordinate vector
    /// @param VecLocD Updated Coordinate
    void setUpdatedCoordinates(VecLocD& coor){coordUpdated_ = coor;};

    /// Updates node coordinate vector
    /// @param int direction @param double updated value
    void updateCoordinate(int dir, double val){coord_(dir) = val;};

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushInverseIncidence(int el) {invIncidence.push_back(el);}

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfElements(){return invIncidence.size();}

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int getInverseIncidenceElement(int i){return invIncidence[i];}
    
    //............................Velocity functions............................
    /// Sets the velocity vector
    /// @param double* velocity vector
    void setVelocity(double *u);

    void setVelocityComponent(int dir, double val){velocity_(dir) = val;};

    void setPreviousVelocityComponent(int dir, double val){previousVelocity_(dir) = val;};

    /// Sets the previous velocity vector
    /// @param double* previous time step velocity vector
    void setPreviousVelocity(double *u);

    /// Sets the intermediate velocity vector
    /// @param double* intermediate velocity vector
    void setIntermediateVelocity(double *u);

    /// Increment the velocity vector
    /// @param int direction @param double increment value
    void incrementVelocity(int dir, double u);

    /// Returns the node velocity vector
    /// @return node velocity vector
    double getVelocity(int dir) {return velocity_(dir);}

    /// Returns the node previous time step velocity vector
    /// @return node previous time step velocity vector
    double getPreviousVelocity(int dir) {return previousVelocity_(dir);}

    /// Returns the node intermediate velocity vector
    /// @return node intermediate velocity vector
    double getIntermediateVelocity(int dir) {return intermediateVelocity_(dir);}

    /// Sets the vorticity at the node
    /// @param double vorticity
    void setVorticity(double div) {vorticity_ = div;}

    /// Returns the nodal vorticity
    /// @return node vorticity
    double getVorticity() {return vorticity_;}

    /// Clears the nodal vorticity
    void clearVorticity() {vorticity_ = 0.;}

    /// Increment the value of the nodal vorticity
    /// @param double increment 
    void incrementVorticity(double val) {if (invIncidence.size() > 0) vorticity_ += val/invIncidence.size();}

    //..........................Acceleration functions..........................
    /// Sets the acceleration vector
    /// @param double* acceleration vector
    void setAcceleration(double *u);

    /// Sets the previous time step acceleration vector
    /// @param double* previous time step acceleration vector
    void setPreviousAcceleration(double *u);

    /// Sets the intermediate acceleration vector
    /// @param double* intermediate acceleration vector
    void setIntermediateAcceleration(double *u);    

    /// Increment the acceleration vector
    /// @param int direction @param double increment value
    void incrementAcceleration(int dir, double u);    

    /// Gets the acceleration vector
    /// @return acceleration vector
    double getAcceleration(int dir) {return acceleration_(dir);}

    /// Gets the previous time step acceleration vector
    /// @return previous time step acceleration vector
    double getPreviousAcceleration(int dir) {return previousAcceleration_(dir);}

    /// Gets the intermediate acceleration vector
    /// @return intermediate acceleration vector
    double getIntermediateAcceleration(int dir) {return intermediateAcceleration_(dir);}

    //............................Pressure functions............................
    /// Sets the nodal pressure
    /// @param double pressure
    void setPressure(double p);

    /// Increments the nodal pressure
    /// @param double pressure increment
    void incrementPressure(double p);

    /// Gets the nodal pressure value
    /// @return nodal pressure value
    double getPressure() {return pressure_;};

    /// Sets the previous time step nodal pressure
    /// @param double previous time step nodal pressure
    void setPreviousPressure(double p);

    /// Gets the previous time step nodal pressure
    /// @return previous time step nodal pressure
    double getPreviousPressure() {return previousPressure_;}

    //.........................Mesh Velocity functions..........................
    /// Sets the node mesh velocity
    /// @param double* mesh velocity
    void setMeshVelocity(double *u);

    /// Sets the previous time step mesh velocity
    /// @param int direction @param double previous time step mesh velocity valu
    void setPreviousMeshVelocity(int dir, double u);

    /// Gets the node mesh velocity
    /// @param int direction @return mesh velocity component
    double getMeshVelocity(int dir) {return meshVelocity_(dir);}

    /// Gets the previous time step mesh velocity
    /// @param int direction @return previous time step mesh velocity component
    double getPreviousMeshVelocity(int dir) {return previousMeshVelocity_(dir);}

    //...........................Constrains functions...........................

    /// Sets all node constrains     
    /// @param int direction 
    /// @param int type: 0 - free, 1 - constrained, 2 - glue zone, 
    /// 3 - fluid-structure interface @param double constrain value
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        velocity_(dir) = value;
        // previousVelocity_(dir) = value;
    };

    /// Gets node constrain type
    /// @return constrain type
    int getConstrains(int dir) {return constrainType[dir];}

    /// Gets node constrain value
    /// @return constrain value
    double getConstrainValue(int dir) {return constrainValue[dir];}

    /// Sets constrains for solving the mesh moving problem
    /// @param int direction 
    /// @param int constrain type: 0 - free, 1 - constrained
    /// @param double constrain value
    void setConstrainsLaplace(int dir, int type, double value){
        constrainTypeLaplace[dir] = type;
        constrainValueLaplace[dir] = value;
    };

    /// Gets constrains of mesh moving problem
    /// @return constrain type
    int getConstrainsLaplace(int dir) {return constrainTypeLaplace[dir];}

    /// Gets constrain value of mesh moving problem
    /// return constrain value
    double getConstrainValueLaplace(int dir){return constrainValueLaplace[dir];}

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------CLEAR VARIABLES--------------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::clearVariables(){
    pressure_ = 0.;          previousPressure_ = 0.;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;   
        intermediateVelocity_(0) = 0.;  
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
        intermediateAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        intermediateVelocity_(1) = 0.;  
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
        intermediateAcceleration_(1) = 0.;
    };
        
    return;
};

template<>
void Node<3>::clearVariables(){

    pressure_ = 0.;          previousPressure_ = 0.;

    if (constrainType[0] != 1){
        velocity_(0) = 0.;   
        previousVelocity_(0) = 0.;  
        intermediateVelocity_(0) = 0.;  
        acceleration_(0) = 0.;
        previousAcceleration_(0) = 0.;
        intermediateAcceleration_(0) = 0.;
    };

    if (constrainType[1] != 1){
        velocity_(1) = 0.;   
        previousVelocity_(1) = 0.;   
        intermediateVelocity_(1) = 0.;
        acceleration_(1) = 0.;
        previousAcceleration_(1) = 0.;
        intermediateAcceleration_(1) = 0.;
    };

    if (constrainType[2] != 1){
        velocity_(2) = 0.;   
        previousVelocity_(2) = 0.;   
        intermediateVelocity_(2) = 0.;
        acceleration_(2) = 0.;
        previousAcceleration_(2) = 0.;
        intermediateAcceleration_(2) = 0.;
    };
        
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setVelocity(double *u){
    //Sets Velocity value
    velocity_(0) = u[0];
    velocity_(1) = u[1]; 
    velocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL VELOCITY VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;
    return;
};

template<>
void Node<3>::incrementVelocity(int dir, double u){
    //All element nodes
    velocity_(dir) += u;

    return;
};


//------------------------------------------------------------------------------
//----------------------INCREMENT NODAL ACCELERATION VALUES---------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementAcceleration(int dir, double u){
    //All element nodes
    acceleration_(dir) += u;
    return;
};

template<>
void Node<3>::incrementAcceleration(int dir, double u){
    //All element nodes
    acceleration_(dir) += u;

    return;
};


//------------------------------------------------------------------------------
//---------------------INCREMENT NODAL COORDINATE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;
    return;
};

template<>
void Node<3>::incrementCoordinate(int dir, double u){
    //All element nodes
    coord_(dir) += u;

    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousCoordinates(int dir, double u){
    //All element nodes
    previousCoord_(dir) = u;
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPressure(double p){
    //Only element nodes 0, 1 and 2
    previousPressure_ = pressure_;  
    pressure_ = p; 
    return;
};

template<>
void Node<3>::setPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    previousPressure_ = pressure_;
    pressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL PRESSURE VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousPressure(double p){
    //Only element nodes 0, 1 and 2  
    previousPressure_ = p; 
    return;
};

template<>
void Node<3>::setPreviousPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    previousPressure_ = p; 
    return;
};

//------------------------------------------------------------------------------
//-----------------------INCREMENT NODAL PRESSURE VALUES------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::incrementPressure(double p){
    //Only element nodes 0, 1 and 2  
    pressure_ += p; 
    return;
};

template<>
void Node<3>::incrementPressure(double p){
    //Only element nodes 0, 1, 2 and 3
    pressure_ += p; 
    return;
};

//------------------------------------------------------------------------------
//---------------------SETS PREVIOUS NODAL VELOCITY VALUES----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousVelocity(double *u){
    //All element nodes
    previousVelocity_(0) = u[0];
    previousVelocity_(1) = u[1]; 
    previousVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//----------------------SETS INTERMEDIATE VELOCITY VALUES-----------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setIntermediateVelocity(double *u){
    //All element nodes
    intermediateVelocity_(0) = u[0];
    intermediateVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setIntermediateVelocity(double *u){
    //All element nodes
    intermediateVelocity_(0) = u[0];
    intermediateVelocity_(1) = u[1]; 
    intermediateVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//--------------------------NODAL ACCELERATION VALUES---------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setAcceleration(double *u){
    //All element nodes
    acceleration_(0) = u[0];
    acceleration_(1) = u[1]; 
    acceleration_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//-------------------SETS PREVIOUS NODAL ACCELERATION VALUES--------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setPreviousAcceleration(double *u){
    //All element nodes
    previousAcceleration_(0) = u[0];
    previousAcceleration_(1) = u[1]; 
    previousAcceleration_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//--------------------SETS INTERMEDIATE ACCELERATION VALUES---------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setIntermediateAcceleration(double *u){
    //All element nodes
    intermediateAcceleration_(0) = u[0];
    intermediateAcceleration_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setIntermediateAcceleration(double *u){
    //All element nodes
    intermediateAcceleration_(0) = u[0];
    intermediateAcceleration_(1) = u[1]; 
    intermediateAcceleration_(2) = u[2]; 
    return;
};


//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setMeshVelocity(double *u){
    //All element nodes
    previousMeshVelocity_(0) = meshVelocity_(0);
    previousMeshVelocity_(1) = meshVelocity_(1); 

    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    return;
};

template<>
void Node<3>::setMeshVelocity(double *u){
    //All element nodes
    meshVelocity_(0) = u[0];
    meshVelocity_(1) = u[1]; 
    meshVelocity_(2) = u[2]; 
    return;
};

//------------------------------------------------------------------------------
//----------------------------NODAL VELOCITY VALUES-----------------------------
//------------------------------------------------------------------------------
template<>
void Node<2>::setPreviousMeshVelocity(int dir, double u){
    //All element nodes
    previousMeshVelocity_(dir) = u;

    return;
};

#endif

