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

#ifndef BOUND_INTEG_QUADRATURE_H
#define BOUND_INTEG_QUADRATURE_H

#include "BoundaryShapeFunction.hpp"

/// Defines the gauss quadrature rule for the boundary integration

template<int DIM>
class BoundaryIntegQuadrature{
public:
    /// Defines the type "PointCoord" which stores the 
    /// integration points coordinates
    typedef ublas::bounded_vector<double, DIM>  PointCoord;

    /// Defines the type "PointWeight" which stores the
    /// integration points weights
    typedef ublas::bounded_vector<double, DIM>      PointWeight;

    /// Defines the numerical integration iterator
    typedef typename PointWeight::iterator              QuadratureListIt;

private:
    //List of integration points coordinates
    PointCoord pointCoord;

    //List of integration points weights
    PointWeight pointWeight;

    double pi = M_PI;

public:
    /// Returns the index of the first integration point
    /// @return First integration point index
    QuadratureListIt begin() {
        return pointWeight.begin();
    }

    /// Returns the index of the last integration point
    /// @return Last integration point index
    QuadratureListIt end() {
        return pointWeight.end();
    }

    /// Returns the integration point adimensional coordinate and weight
    /// @return Adimensional coordinate and weight of the integration point
    std::pair<PointCoord,PointWeight> GaussQuadrature(){

        double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
        int mga;
        double nga = DIM;

        mga=(nga+1.)/2.;
        xmga=0.0;
        xlga=1.0;
        
        for (int iga=1; iga<=mga; iga++) {
            zga = cos(pi*(double(iga)-0.25)/(double(nga)+0.5));
        g1:
            p1ga = 1.0;
            p2ga = 0.0;
            for (int jga=1; jga <= nga; jga++) {
                p3ga = p2ga;
                p2ga = p1ga;
                p1ga = ((2.0*double(jga)-1.0)*zga*p2ga-     \
                        (double(jga)-1.0)*p3ga)/(double(jga));
            };
         
            ppga = nga*(zga*p1ga-p2ga)/(zga*zga-1.0);
            z1ga = zga;
            zga = z1ga-p1ga/ppga;
            
            if (fabs(zga-z1ga) > 1.0e-15) goto g1;
            
            pointCoord(iga-1) = xmga-xlga*zga;
            pointCoord(nga-iga) = xmga + xlga*zga;
            pointWeight(iga-1) = 2.0*xlga/((1.0-zga*zga)*ppga*ppga);
            pointWeight(nga-iga) = pointWeight(iga-1);             
        };

        return std::make_pair(pointCoord,pointWeight);
    }; 


};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

#endif
