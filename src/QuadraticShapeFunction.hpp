//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------QUADRATIC SHAPE FUNCTION----------------------------
//------------------------------------------------------------------------------

#ifndef QUADSHAPEFUNCTION_H
#define QUADSHAPEFUNCTION_H

using namespace boost::numeric;

/// Defines the quadratic shape functions and its derivatives

template<int DIM>

class QuadShapeFunction {
public:

    /// Defines type "Coords" to allocate the integration points coordinates
    typedef ublas::bounded_vector<double, DIM>     Coords;

    //Number of interpolation notes    
    static const int numIntpNodes = 4*DIM-2;
    
    /// Defines the type "Values" which stores the shape function values
    typedef ublas::bounded_vector<double, numIntpNodes>           Values;

    /// Defines the type "ValuesDeriv" which stores the shape function
    /// derivatives
    typedef ublas::bounded_matrix<double, DIM, numIntpNodes>      ValueDeriv;

    /// Defines the type "ValuesDDeriv" which stores 
    /// the shape function second derivatives
    typedef ublas::bounded_vector<double, numIntpNodes>           ValueDDphi;
    typedef ublas::bounded_matrix<ValueDDphi, DIM, DIM>           ValueDDeriv;

public:
    
    /// Evaluates the shape function value
    /// @param Coords Adimensional coordinates 
    /// @param Values Shape function values
    void evaluate(const Coords& xi, Values& phi) const;
    
    /// Evaluates the values of the shape funtion derivatives
    /// @param Coords Adimensional coordinates 
    /// @param ValueDeriv Shape function derivatives values
    void evaluateGradient(const Coords& xi, ValueDeriv& dphi) const;   

    /// Evaluates the values of the shape funtion second derivatives    
    /// @param Coords Adimensional coordinates 
    /// @param ValueDeriv Shape function second derivatives values
    void evaluateHessian(const Coords& xi, ValueDDeriv& ddphi) const;   
};


//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void QuadShapeFunction<2>::evaluate(const Coords& xi, Values& phi) const {
    
     const double xsi1 = xi(0);
     const double xsi2 = xi(1);
     const double xsi3 = 1. - xsi1 - xsi2;

     phi(0) = xsi3 * (2.0 * xsi3 - 1.0);
     phi(1) = xsi1 * (2.0 * xsi1 - 1.0);
     phi(2) = xsi2 * (2.0 * xsi2 - 1.0);
     phi(3) = 4.0 * xsi3 * xsi1;
     phi(4) = 4.0 * xsi1 * xsi2;
     phi(5) = 4.0 * xsi2 * xsi3;

     // element conectivity
     //     2
     //     54
     //     031
     
     return;
}

template<>
void QuadShapeFunction<3>::evaluate(const Coords& xi, Values& phi) const {
    
     const double xsi1 = xi(0);
     const double xsi2 = xi(1);
     const double xsi3 = xi(2);

     phi(0) = 2.0 * (xsi1 - 0.50) * xsi1;
     phi(1) = 2.0 * (xsi2 - 0.50) * xsi2;
     phi(2) = 2.0 * (xsi3 - 0.50) * xsi3;
     phi(3) = (2.0 - 2.0 * xsi3 - 2.0 * xsi2 - 2.0 * xsi1 - 1.0)  \
         * (1.0 - xsi1 - xsi2 - xsi3);
     phi(4) = 4.0 * xsi1 * xsi2;
     phi(5) = 4.0 * xsi1 * xsi3;
     phi(6) = 4.0 * xsi1 * (1.0 - xsi1 - xsi2 - xsi3);
     phi(7) = 4.0 * xsi2 * xsi3;
     phi(8) = 4.0 * xsi3 * (1.0 - xsi1 - xsi2 - xsi3);
     phi(9) = 4.0 * xsi2 * (1.0 - xsi1 - xsi2 - xsi3);

     // element conectivity
     //layer 1
     //     3
     //     8 9
     //     2 7 1
     //layer 2
     //     6
     //     54
     //layer 3
     //     0
     
     return;
}

//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2>::evaluateGradient(const Coords& xi, \
                                            ValueDeriv& dphi) const {

     const double xsi1 = xi(0);
     const double xsi2 = xi(1);
     const double xsi3 = 1. - xsi1 - xsi2;

     dphi(0,1) = 4. * xsi1 - 1.;
     dphi(1,1) = 0.;

     dphi(0,2) = 0.;
     dphi(1,2) = 4. * xsi2 - 1.;
 
     dphi(0,0) = -4. * xsi3 + 1.;
     dphi(1,0) = -4. * xsi3 + 1.;

     dphi(0,4) = 4. * xsi2;
     dphi(1,4) = 4. * xsi1;

     dphi(0,5) = -4. * xsi2;
     dphi(1,5) = 4. * (xsi3 - xsi2);

     dphi(0,3) = 4. * (xsi3 - xsi1);
     dphi(1,3) = -4. * xsi1;

     // element conectivity
     //     2
     //     54
     //     031
     
    return;
}

template<>
void QuadShapeFunction<3>::evaluateGradient(const Coords& xi, \
                                            ValueDeriv& dphi) const {

     const double xsi1 = xi(0);
     const double xsi2 = xi(1);
     const double xsi3 = xi(2);
     
     dphi(0,0) = 4. * xsi1 - 1.;
     dphi(1,0) = 0.;
     dphi(2,0) = 0.;

     dphi(0,1) = 0.;
     dphi(1,1) = 4. * xsi2 - 1.;
     dphi(2,1) = 0.;
     
     dphi(0,2) = 0.;
     dphi(1,2) = 0.;
     dphi(2,2) = 4. * xsi3 - 1.;

     dphi(0,3) = 4. * (xsi1 + xsi2 + xsi3) - 3.;
     dphi(1,3) = 4. * (xsi1 + xsi2 + xsi3) - 3.;
     dphi(2,3) = 4. * (xsi1 + xsi2 + xsi3) - 3.;

     dphi(0,4) = 4. * xsi2;
     dphi(1,4) = 4. * xsi1;
     dphi(2,4) = 0.;

     dphi(0,5) = 4. * xsi3;
     dphi(1,5) = 0.;
     dphi(2,5) = 4. * xsi1;

     dphi(0,6) = 4. * (1. - 2. * xsi1 - xsi2 - xsi3);
     dphi(1,6) = -4. * xsi1;
     dphi(2,6) = -4. * xsi1;

     dphi(0,7) = 0.;
     dphi(1,7) = 4. * xsi3;
     dphi(2,7) = 4. * xsi2;

     dphi(1,8) = -4. * xsi3;
     dphi(0,8) = -4. * xsi3;
     dphi(2,8) = 4. * (1. - 2. * xsi3 - xsi2 - xsi1);

     dphi(0,9) = -4. * xsi2;
     dphi(1,9) = 4. * (1. - 2. * xsi2 - xsi1 - xsi3);
     dphi(2,9) = -4. * xsi2;

     // element conectivity
     //layer 1
     //     3
     //     8 9
     //     2 7 1
     //layer 2
     //     6
     //     54
     //layer 3
     //     0
     
    return;
}

//------------------------------------------------------------------------------
//----------------COMPUTE SHAPE FUNCTION SECOND DERIVATIVE VALUE----------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2>::evaluateHessian(const Coords& xi, \
                                           ValueDDeriv& ddphi) const {

     ddphi(0,0)(0) = 4.;
     ddphi(0,1)(0) = 4.;
     ddphi(1,0)(0) = 4.;
     ddphi(1,1)(0) = 4.;

     ddphi(0,0)(1) = 4.;
     ddphi(0,1)(1) = 0.;
     ddphi(1,0)(1) = 0.;
     ddphi(1,1)(1) = 0.;

     ddphi(0,0)(2) = 0.;
     ddphi(0,1)(2) = 0.;
     ddphi(1,0)(2) = 0.;
     ddphi(1,1)(2) = 4.;

     ddphi(0,0)(3) = -8.;
     ddphi(0,1)(3) = -4.;
     ddphi(1,0)(3) = -4.;
     ddphi(1,1)(3) = 0.;

     ddphi(0,0)(4) = 0.;
     ddphi(0,1)(4) = 4.;
     ddphi(1,0)(4) = 4.;
     ddphi(1,1)(4) = 0.;

     ddphi(0,0)(5) = 0.;
     ddphi(0,1)(5) = -4.;
     ddphi(1,0)(5) = -4.;
     ddphi(1,1)(5) = -8.;

     // element conectivity
     //     2
     //     54
     //     031
     
    return;
}

template<>
void QuadShapeFunction<3>::evaluateHessian(const Coords& xi, \
                                           ValueDDeriv& ddphi) const {

     ddphi(0,0)(0) = 4.;
     ddphi(0,1)(0) = 0.;
     ddphi(0,2)(0) = 0.;
     ddphi(1,0)(0) = 0.;
     ddphi(1,1)(0) = 0.;
     ddphi(1,2)(0) = 0.;
     ddphi(2,0)(0) = 0.;
     ddphi(2,1)(0) = 0.;
     ddphi(2,2)(0) = 0.;

     ddphi(0,0)(1) = 0.;
     ddphi(0,1)(1) = 0.;
     ddphi(0,2)(1) = 0.;
     ddphi(1,0)(1) = 0.;
     ddphi(1,1)(1) = 4.;
     ddphi(1,2)(1) = 0.;
     ddphi(2,0)(1) = 0.;
     ddphi(2,1)(1) = 0.;
     ddphi(2,2)(1) = 0.;

     ddphi(0,0)(2) = 0.;
     ddphi(0,1)(2) = 0.;
     ddphi(0,2)(2) = 0.;
     ddphi(1,0)(2) = 0.;
     ddphi(1,1)(2) = 0.;
     ddphi(1,2)(2) = 0.;
     ddphi(2,0)(2) = 0.;
     ddphi(2,1)(2) = 0.;
     ddphi(2,2)(2) = 4.;

     ddphi(0,0)(3) = 4.;
     ddphi(0,1)(3) = 4.;
     ddphi(0,2)(3) = 4.;
     ddphi(1,0)(3) = 4.;
     ddphi(1,1)(3) = 4.;
     ddphi(1,2)(3) = 4.;
     ddphi(2,0)(3) = 4.;
     ddphi(2,1)(3) = 4.;
     ddphi(2,2)(3) = 4.;

     ddphi(0,0)(4) = 0.;
     ddphi(0,1)(4) = 4.;
     ddphi(0,2)(4) = 0.;
     ddphi(1,0)(4) = 4.;
     ddphi(1,1)(4) = 0.;
     ddphi(1,2)(4) = 0.;
     ddphi(2,0)(4) = 0.;
     ddphi(2,1)(4) = 0.;
     ddphi(2,2)(4) = 0.;

     ddphi(0,0)(5) = 0.;
     ddphi(0,1)(5) = 0.;
     ddphi(0,2)(5) = 4.;
     ddphi(1,0)(5) = 0.;
     ddphi(1,1)(5) = 0.;
     ddphi(1,2)(5) = 0.;
     ddphi(2,0)(5) = 4.;
     ddphi(2,1)(5) = 0.;
     ddphi(2,2)(5) = 0.;

     ddphi(0,0)(6) = -8.;
     ddphi(0,1)(6) = -4.;
     ddphi(0,2)(6) = -4.;
     ddphi(1,0)(6) = -4.;
     ddphi(1,1)(6) = 0.;
     ddphi(1,2)(6) = 0.;
     ddphi(2,0)(6) = -4.;
     ddphi(2,1)(6) = 0.;
     ddphi(2,2)(6) = 0.;

     ddphi(0,0)(7) = 0.;
     ddphi(0,1)(7) = 0.;
     ddphi(0,2)(7) = 0.;
     ddphi(1,0)(7) = 0.;
     ddphi(1,1)(7) = 0.;
     ddphi(1,2)(7) = 4.;
     ddphi(2,0)(7) = 0.;
     ddphi(2,1)(7) = 4.;
     ddphi(2,2)(7) = 0.;

     ddphi(0,0)(8) = 0.;
     ddphi(0,1)(8) = 0.;
     ddphi(0,2)(8) = -4.;
     ddphi(1,0)(8) = 0.;
     ddphi(1,1)(8) = 0.;
     ddphi(1,2)(8) = -4.;
     ddphi(2,0)(8) = -4.;
     ddphi(2,1)(8) = -4.;
     ddphi(2,2)(8) = -8.;

     ddphi(0,0)(9) = 0.;
     ddphi(0,1)(9) = -4.;
     ddphi(0,2)(9) = 0.;
     ddphi(1,0)(9) = -4.;
     ddphi(1,1)(9) = -8.;
     ddphi(1,2)(9) = -4.;
     ddphi(2,0)(9) = 0.;
     ddphi(2,1)(9) = -4.;
     ddphi(2,2)(9) = 0.;

//      // element conectivity
//      //layer 1
//      //     3
//      //     8 9
//      //     2 7 1
//      //layer 2
//      //     6
//      //     54
//      //layer 3
//      //     0
     
    return;
}



#endif
