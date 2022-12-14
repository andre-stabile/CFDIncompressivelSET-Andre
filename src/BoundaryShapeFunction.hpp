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

#ifndef BOUND_SHAPEFUNCTION_H
#define BOUND_SHAPEFUNCTION_H

using namespace boost::numeric;

/// Defines the fluid boundary shape functions

template<int DIM>
class BoundShapeFunction {
public:
   
    /// Defines the type "Values" which stores the shape function values
    typedef ublas::bounded_vector<double, 3*(DIM-1)>           Values;

    /// Defines the type "ValuesDeriv" which stores the shape function
    /// derivatives
    typedef ublas::bounded_vector<double, 3*(DIM-1)>           ValueDeriv;

public:
    
    /// Evaluates the shape function value
    /// @param double Adimensional coordinate
    void evaluate(double Xsi);

    /// Returns the shape function value
    /// @param double Adimensional coordinates 
    /// @return boundary shape function value
    Values getShapeFunction(double Xsi){
        evaluate(Xsi);
        return phi_;
    };

    /// Returns the shape function derivative value
    /// @param double Adimensional coordinates 
    /// @return boundary shape function derivative value
    Values getShapeFunctionDerivative(double Xsi){
        evaluate(Xsi);
        return dphi_;
    };

private:
    Values     phi_;
    ValueDeriv dphi_;
    
};


//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void BoundShapeFunction<2>::evaluate(double Xsi){
    
    double aux;
    double Nnos = 3;
    double AdimCoord[3];
    AdimCoord[0] = -1.;
    AdimCoord[1] =  0.;
    AdimCoord[2] =  1.;

    
    for (int j=0; j<Nnos; j++) {

        phi_(j) = 1.0;
        dphi_(j) = 0.0;

        for (int i=0; i<Nnos; i++) {

            aux=1.0;

            if(i != j){
                
                phi_(j) = phi_(j) * (Xsi - AdimCoord[i]) / 
                    (AdimCoord[j] - AdimCoord[i]);
                
                for (int k=0; k<Nnos; k++) {
                    if ((i != k) && (j != k)) aux = aux*(Xsi-AdimCoord[k]);
                };
                dphi_(j) += aux;
            };
        };
    };

    for (int i=0; i<Nnos; i++) {
        for (int k=0; k<Nnos; k++) {
            if (i != k) dphi_(i) = dphi_(i)/(AdimCoord[i]-AdimCoord[k]);
        };
    };    
    
    return;
}




#endif
