//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef INTEG_QUADRATURE_H
#define INTEG_QUADRATURE_H

#include "QuadraticShapeFunction.hpp"

/// Defines the domain integration Hammer quadrature

template<int DIM>
class IntegQuadrature{
public:
    /// Defines the type "PointCoord" which stores the 
    /// integration points coordinates
    typedef ublas::bounded_matrix<double, 4*DIM-1,DIM>  PointCoord;

    /// Defines the type "PointWeight" which stores the
    /// integration points weights
    typedef ublas::bounded_vector<double, 4*DIM-1>      PointWeight;

    /// Integration point logical vector
    typedef ublas::bounded_vector<bool, 4*DIM-1>        PointLogical;

    /// integration points int variable
    typedef ublas::bounded_vector<int, 4*DIM-1>         PointIntVar;

    /// Defines the numerical integration iterator
    typedef typename PointWeight::iterator              QuadratureListIt;

    /// Defines vector of nodal values
    typedef ublas::bounded_vector<double, 4*DIM-2>      NodalValuesQuad;

public:
    /// Returns the index of the first integration point
    /// @return first integration point index
    QuadratureListIt begin() {
        return pointWeight.begin();
    }

    /// Returns the index of the last integration point
    /// @return last integration point index
    QuadratureListIt end() {
        return pointWeight.end();
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j); 

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i);
  
    /// Interpolate quadratic variables
    /// @param NodalValuesQuad element variable nodal values
    /// @param Integration point index
    /// @return Interpolated variable value
    double interpolateQuadraticVariable(NodalValuesQuad nValues, int point);


private:
    //List of integration points coordinates
    PointCoord pointCoord;

    //List of integration points weights
    PointWeight pointWeight;

    //Defines shape functions
    QuadShapeFunction<DIM> shapeQuad;

    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::PointList(int i, int j){

    pointCoord(0,0) = 1. / 3.;
    pointCoord(0,1) = 1. / 3.;
        
    pointCoord(1,0) = (9. + 2. * sqrt(15.)) / 21.;
    pointCoord(1,1) = (6. - sqrt(15.)) / 21.;
      
    pointCoord(2,0) = (6. - sqrt(15.)) / 21.;
    pointCoord(2,1) = (9. + 2. * sqrt(15.)) / 21.;
      
    pointCoord(3,0) = (6. - sqrt(15.)) / 21.;
    pointCoord(3,1) = (6. - sqrt(15.)) / 21.;
      
    pointCoord(4,0) = (6. + sqrt(15.)) / 21.;
    pointCoord(4,1) = (6. + sqrt(15.)) / 21.;
      
    pointCoord(5,0) = (9. - 2. * sqrt(15.)) / 21.;
    pointCoord(5,1) = (6. + sqrt(15.)) / 21.;
      
    pointCoord(6,0) = (6. + sqrt(15.)) / 21.;
    pointCoord(6,1) = (9. - 2. * sqrt(15.)) / 21.;

 
    // pointCoord(0,0) = 1. / 3.;
    // pointCoord(0,1) = 1. / 3.;
        
    // pointCoord(1,0) = 0.797426985353087;
    // pointCoord(1,1) = 0.101286507323456;
      
    // pointCoord(2,0) = 0.101286507323456;
    // pointCoord(2,1) = 0.797426985353087;
      
    // pointCoord(3,0) = 0.101286507323456;
    // pointCoord(3,1) = 0.101286507323456;
      
    // pointCoord(4,0) = 0.470142064105115;
    // pointCoord(4,1) = 0.470142064105115;
      
    // pointCoord(5,0) = 0.05971587178977;
    // pointCoord(5,1) = 0.470142064105115;
      
    // pointCoord(6,0) = 0.470142064105115;
    // pointCoord(6,1) = 0.05971587178977;

    return pointCoord(i,j);
};

template<>
double IntegQuadrature<3>::PointList(int i, int j){
    
    const double a = (1. + sqrt(5. / 14.)) / 4.;
    const double b = (1. - sqrt(5. / 14.)) / 4.;

    pointCoord(0,0) = 1. / 4.;
    pointCoord(0,1) = 1. / 4.;
    pointCoord(0,2) = 1. / 4.;

    pointCoord(1,0) = 11. / 14.;
    pointCoord(1,1) = 1. / 14.;
    pointCoord(1,2) = 1. / 14.;

    pointCoord(2,0) = 1. / 14.;
    pointCoord(2,1) = 11. / 14.;
    pointCoord(2,2) = 1. / 14.;

    pointCoord(3,0) = 1. / 14.;
    pointCoord(3,1) = 1. / 14.;
    pointCoord(3,2) = 11. / 14.;

    pointCoord(4,0) = 1. / 14.;
    pointCoord(4,1) = 1. / 14.;
    pointCoord(4,2) = 1. / 14.;

    pointCoord(5,0) = a;
    pointCoord(5,1) = a;
    pointCoord(5,2) = b;

    pointCoord(6,0) = a;
    pointCoord(6,1) = b;
    pointCoord(6,2) = a;

    pointCoord(7,0) = a;
    pointCoord(7,1) = b;
    pointCoord(7,2) = b;

    pointCoord(8,0) = b;
    pointCoord(8,1) = a;
    pointCoord(8,2) = a;

    pointCoord(9,0) = b;
    pointCoord(9,1) = a;
    pointCoord(9,2) = b;

    pointCoord(10,0) = b;
    pointCoord(10,1) = b;
    pointCoord(10,2) = a;

    return pointCoord(i,j);
}

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::WeightList(int i){
    
    pointWeight(0) = 0.11250;
    pointWeight(1) = (155. - sqrt(15.)) / 2400.;
    pointWeight(2) = (155. - sqrt(15.)) / 2400.;
    pointWeight(3) = (155. - sqrt(15.)) / 2400.;
    pointWeight(4) = (155. + sqrt(15.)) / 2400.;
    pointWeight(5) = (155. + sqrt(15.)) / 2400.;
    pointWeight(6) = (155. + sqrt(15.)) / 2400.; 

    // pointWeight(0) = 0.11250;
    // pointWeight(1) = 0.125939180544827 / 2.;
    // pointWeight(2) = 0.125939180544827 / 2.;
    // pointWeight(3) = 0.125939180544827 / 2.;
    // pointWeight(4) = 0.132394152788506 / 2.;
    // pointWeight(5) = 0.132394152788506 / 2.;
    // pointWeight(6) = 0.132394152788506 / 2.; 

    return pointWeight(i);
};

template<>
double IntegQuadrature<3>::WeightList(int i){
    
    pointWeight(0) = -74. / 5625.;
    pointWeight(1) = 343. / 45000.;
    pointWeight(2) = 343. / 45000.;
    pointWeight(3) = 343. / 45000.;
    pointWeight(4) = 343. / 45000.;
    pointWeight(5) = 56. / 2250.;
    pointWeight(6) = 56. / 2250.; 
    pointWeight(7) = 56. / 2250.;
    pointWeight(8) = 56. / 2250.; 
    pointWeight(9) = 56. / 2250.;
    pointWeight(10) = 56. / 2250.; 

    return pointWeight(i);
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 2> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 6; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};

template<>
double IntegQuadrature<3>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 3> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    xsi(2) = PointList(point,2);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 10; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};


#endif
