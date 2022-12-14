//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------------FLUID BOUNDARY--------------------------------
//------------------------------------------------------------------------------

#ifndef BOUNDARY_H
#define BOUNDARY_H

/// Defines the fluid boundary object and its properties

template<int DIM>
class Boundary{

public: 
    /// Defines a type to store the boundary element connectivity vector
    typedef ublas::bounded_vector<int, 3*(DIM-1)>         BoundConnect;

    /// Boundary element constructor
    /// @param BoundConnect boundary element connectivity
    /// @param int boundary element index
    /// @param int boundary element constrain type 
    /// @param int boundary element constrain value @see Node::setConstrains()
    Boundary(BoundConnect& connec, int index,  \
             int *constrain, double *values, int gr){
        connectB_ = connec;
        index_ = index;
        group_ = gr;

        constrainType[0] = constrain[0];
        constrainType[1] = constrain[1];
        constrainType[2] = constrain[2];
        
        constrainValue[0] = values[0];
        constrainValue[1] = values[1];
        constrainValue[2] = values[2];

        element_ = 1.e10;
        elementSide_ = 1000;
    };     

    /// Returns the boundary element constrain component type
    /// @param int direction @return boundary element constrain component type
    int getConstrain(int dir){return constrainType[dir];}

    /// Returns the boundary element constrain component value
    /// @param int direction @return boundary element constrain component value
    double getConstrainValue(int dir){return constrainValue[dir];}

    /// Returns the boundary element connectivity
    /// @return boundary element connectivity
    BoundConnect getBoundaryConnectivity(){return connectB_;}

    /// Sets the boundary element group
    /// @param int boundary element group
    void setBoundaryGroup(int gr){group_ = gr;}

    /// Gets the boundary element group
    /// @return boundary element group
    int getBoundaryGroup(){return group_;};

    /// Sets the fluid element correspondence
    /// @param int fluid element correspondence
    void setElement(int el){element_ = el;}

    /// Sets the fluid element correspondence side at the boundary
    /// @param int fluid element correspondence side
    void setElementSide(int el){elementSide_ = el;}

    /// Gets the fluid element correspondence
    /// @return fluid element correspondence
    int getElement(){return element_;}

    /// Gets the fluid element correspondence side at the boundary
    /// @return fluid element correspondence side
    int getElementSide(){return elementSide_;}

private:
    BoundConnect connectB_;         //Boundary element connectivity
    int          index_;            //Boundary element index
    int          constrainType[3];  //Element type of constrain
    double       constrainValue[3]; //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group
};






















#endif
