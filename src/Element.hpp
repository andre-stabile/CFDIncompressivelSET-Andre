//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//------------------------------------------------------------------------------

#ifndef ELEMENT_H
#define ELEMENT_H

#include "Node.hpp"
#include "BoundaryIntegrationQuadrature.hpp"
#include "IntegrationQuadrature.hpp"

/// Defines the fluid element object and all the element information
template<int DIM>
class Element{
 
public:
    /// Defines the class Node locally
    typedef Node<DIM>                                           Nodes;

    /// Defines a type to store the element mesh connectivity
    typedef ublas::bounded_vector<int, 4*DIM-2>                 Connectivity;

    /// Defines a blas-type vector with dimension = DIM
    typedef ublas::bounded_vector<double, DIM>                  DimVector;

    /// Defines a blas-type matrix with dimension = DIM x DIM
    typedef ublas::bounded_matrix<double, DIM, DIM>             DimMatrix;

    /// Defines the vector which contains the element nodal coordinates
    typedef ublas::bounded_matrix<double, 4*DIM-2, DIM>         LocalNodes;

    /// Defines the local vector type with dimension 18 for DIM=2 
    /// and 40 for DIM=3
    typedef ublas::bounded_vector<double, 22*DIM-26>            LocalVector;

    /// Defines the local matrix type with dimension 18x18
    /// for DIM=2 and 40x40 for DIM=3
    typedef ublas::bounded_matrix<double, 22*DIM-26, 22*DIM-26> LocalMatrix;

    /// Defines the normal integration quadrature rule class locally
    typedef IntegQuadrature<DIM>                                NormalQuad;

    /// Defines the boundary integration quadrature rule class locally
    typedef BoundaryIntegQuadrature<DIM>                        BoundaryQuad;
    
    /// Define the type VecLocD from class Node locally
    typedef typename Nodes::VecLocD                             VecLoc;

private:
    QuadShapeFunction<DIM> shapeQuad; //Quadratic shape function
    BoundShapeFunction<DIM>shapeBound;//Boundary shape function
    NormalQuad             nQuad;     //Integration quadrature
    std::vector<Nodes *>   nodes_;    //Velocity nodes
    Connectivity  connect_;           //Velocity mesh connectivity 
    int           index_;             //Element index
    LocalNodes    localNodesBoundary_;//Nodal coordinates - velocity
    double        djac_;              //Jacobian determinant
    double        weight_;            //Integration point weight
    DimMatrix     ainv_;              //Inverse of Jacobian transform
    double        visc_;              //Fluid dynamic viscosity
    double        dens_;              //Fluid Density
    double        dTime_;             //Time Step
    double        timeScheme_;        //Time Integration scheme
    LocalMatrix   jacobianNRMatrix;   //Newton's method jacobian
    double        x_, y_;
    double        u_, v_, w_, p_;     //Interpolated velocity and pressure
    double        uPrev_, vPrev_, wPrev_, pPrev_;
    double        uInt_, vInt_, axInt_, ayInt_;
    double        ax_, ay_, az_;      //Interpolated acceleration
    double        umesh_, vmesh_, wmesh_; //Interpolated mesh velocity
    double        du_dx, du_dy, dv_dx, dv_dy; //Interpolated fluid spatial derivatives
    double        duprev_dx, duprev_dy, dvprev_dx, dvprev_dy;
    double        dax_dx, dax_dy, day_dx, day_dy;
    double        du_dxx, du_dyy, du_dxy, dv_dxx, dv_dyy, dv_dxy;
    double        dp_dx, dp_dy, dp_dz;//Interpolated pressure spatial derivative
    double        dumesh_dx, dumesh_dy, dumesh_dz,//Interpolated mesh 
                  dvmesh_dx, dvmesh_dy, dvmesh_dz,//velocity derivatives
                  dwmesh_dx, dwmesh_dy, dwmesh_dz;
    double        dpPrev_dx, dpPrev_dy;
    double        tSUPG_;             //SUPG parameter
    double        tPSPG_;             //PSPG parameter
    double        tLSIC_;             //LSIC parameter
    double        tSUGN1_;
    double        tSUGN2_;
    double        tSUGN3_;
    double        hRGN_;
    DimVector     fieldForce_;        //Element field force
    LocalVector   externalForces;     //Field forces integrated over
                                      //the element domain
    LocalVector   rhsVector;          //RHS vector of Newton's method
    int           sideBoundary_;
    double        meshMovingParameter;
    double        pressureDragForce;
    double        pressureLiftForce;
    double        frictionDragForce;
    double        frictionLiftForce;
    double        dragForce;
    double        liftForce;
    double        pitchingMoment;
    double        perimeter;

    //Second derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDDeriv ddphi_dx; 
    //First derivatives of velocity shape functions
    typename QuadShapeFunction<DIM>::ValueDeriv  dphi_dx, dphiL_dx;
    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::Values     phib_;     
    //Values of velocity shape functins
    typename BoundShapeFunction<DIM>::ValueDeriv dphib_;     

public:
    /// fluid element constructor
    /// @param int element index @param Connectivity element connectivity
    /// @param vector<Nodes> 
    Element(int index, Connectivity& connect, std::vector<Nodes *> &nodes){
        index_ = index;
        connect_ = connect;
        nodes_ = nodes;

        djac_ = 0.;             weight_ = 0.;             ainv_.clear();
        visc_ = 0.;             dens_ = 0.;               dTime_ = 0.; 
        timeScheme_ = 0.;   
        jacobianNRMatrix.clear();
        umesh_ = 0.;      vmesh_ = 0.;
        fieldForce_.clear();     externalForces.clear();     rhsVector.clear();
      
        sideBoundary_ = 0;

        DimVector xsi;        
        xsi.clear();        
        getJacobianMatrix(xsi);
    };

    //........................Element basic information.........................
    /// Clear all element variables
    void clearVariables();

    /// Sets the element connectivity
    /// @param Connectivity element connectivity
    void setConnectivity(Connectivity& connect){connect_ = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    Connectivity getConnectivity(){return connect_;};

    /// Sets the element density
    /// @param double element density
    void setDensity(double& dens){dens_ = dens;}

    /// Sets the element viscosity
    /// @param double element viscosity
    void setViscosity(double& visc){visc_ = visc;}

    /// Sets the time step size
    /// @param double time step size
    void setTimeStep(double& dt){dTime_ = dt;}

    /// Sets the time integration scheme
    /// @param double time integration scheme: 0.0 - Explicit forward Euler;
    /// 1.0 - Implicit backward Euler;
    /// 0.5 - Implicit Trapezoidal Rule.
    void setTimeIntegrationScheme(double& b){timeScheme_ = b;}

    /// Sets the body forces
    /// @param double* body forces
    void setFieldForce(double* ff);

    /// Returns the element index
    int getIndex(){return index_;};

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    void getJacobianMatrix(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    void getSpatialDerivatives(ublas::bounded_vector<double, DIM>& xsi);

    /// Compute and stores the interpolated velocities, mesh velocities, 
    /// previous mesh velocity, acceleration, among others
    void getVelAndDerivatives();

    /// Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getParameterSUPG();
    double getPSPG(){return tPSPG_;};

    /// Compute the vorticity field
    void computeVorticity();

    /// Gets the element jacobian determinant
    /// @return element jacobinan determinant
    double getJacobian(){return djac_;};

    /// Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces();

    /// Gets the element pressure drag force
    double getPressureDragForce(){return pressureDragForce;};

    /// Gets the element pressure lift force
    double getPressureLiftForce(){return pressureLiftForce;};

    /// Gets the element friction drag force
    double getFrictionDragForce(){return frictionDragForce;};

    /// Gets the element friction lift force
    double getFrictionLiftForce(){return frictionLiftForce;};

    /// Gets the element drag force
    double getDragForce(){return dragForce;};

    /// Gets the element lift force
    double getLiftForce(){return liftForce;};

    /// Gets the element Pitching Moment
    double getPitchingMoment(){return pitchingMoment;};

    /// Gets the element boundary perimeter
    double getPerimeter(){return perimeter;};

    /// Gets the spatial jacobian matrix
    /// @param bounded_vector integration point coordinates
    /// @return spatial jacobian matrix
    DimMatrix getJacobianMatrixValues(ublas::bounded_vector<double, DIM>& xsi){
        getJacobianMatrix(xsi);
        return ainv_;};

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){sideBoundary_ = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int getElemSideInBoundary(){return sideBoundary_;};

    /// Sets the mesh moving weighting parameter for solving the Laplace problem
    /// @param double parameter value
    void setMeshMovingParameter(double value) {meshMovingParameter = value;};

    /// Gets the mesh moving weighting parameter
    /// @return mesh moving weighting parameter
    double getMeshMovingParameter(){return meshMovingParameter;};

    //.......................Element vectors and matrices.......................
    /// Compute and store the element matrix for the incompressible flow problem
    /// @param int integration point index
    void getElemMatrix(int index);

    /// Compute and store the element matrix for the Laplace/Poisson problem
    void getElemLaplMatrix();

    /// Sets the boundary conditions for the incompressible flow problem
    void setBoundaryConditions();

    /// Sets the boundary conditions for the Laplace/Poisson problem
    void setBoundaryConditionsLaplace();

    ///Compute and store the residual vector for the incompressible flow problem
    /// @param int integration point index
    void getResidualVector(int index);

    /// Compute and store the residual vector for the Laplace/Poisson problem
    void getResidualVectorLaplace();

    /// Gets the Newton-Raphson's jacobian matrix
    /// @return Newton-Raphson's jacobian matrix
    LocalMatrix getJacNRMatrix(){return jacobianNRMatrix;};

    /// Gets the residual vector
    /// @return residual vecot
    LocalVector getRhsVector(){return rhsVector;};

    /// Compute and store the nodal gradient value (for potential problems)
    void computeNodalGradient();

    //...............................Problem type...............................
    /// Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes();

    /// Compute the Steady Laplace problem matrices and vectors 
    /// (usually for the mesh moving step)
    void getSteadyLaplace();

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------SET ELEMENT FIELD FORCE---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setFieldForce(double *ff) {
    fieldForce_(0) = ff[0];
    fieldForce_(1) = ff[1];
    
    return;
};

//------------------------------------------------------------------------------
//---------------------------CLEAR ELEMENT VARIABLES----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::clearVariables(){

    djac_ = 0.;               weight_ = 0.;              ainv_.clear();
    
    jacobianNRMatrix.clear();
    
    u_ = 0.;          v_ = 0.;          w_ = 0;          p_ = 0.; 
    umesh_ = 0.;      vmesh_ = 0.;      wmesh_ = 0.; 
    du_dx = 0.;       du_dy = 0.;       
    dv_dx = 0.;       dv_dy = 0.;       
    tSUPG_ = 0.;      tPSPG_ = .0;      tLSIC_ = 0.;
    
    externalForces.clear();   rhsVector.clear();

    ddphi_dx.clear();        dphi_dx.clear();            phi_.clear();      

    return;
}; 

//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getJacobianMatrix(ublas::bounded_vector<double,2>& xsi) {

    //Computes the spatial Jacobian matrix and its inverse

    typename QuadShapeFunction<2>::ValueDeriv dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    // for (int i=0; i<6; i++){
    //     dx_dxsi1 += localNodes_(i,0) * dphi(0,i);
    //     dx_dxsi2 += localNodes_(i,0) * dphi(1,i);
    //     dy_dxsi1 += localNodes_(i,1) * dphi(0,i);
    //     dy_dxsi2 += localNodes_(i,1) * dphi(1,i);        
    // };

    for (int i=0; i<6; i++){
        dx_dxsi1 += nodes_[connect_(i)] -> getCoordinateValue(0) * dphi(0,i);
        dx_dxsi2 += nodes_[connect_(i)] -> getCoordinateValue(0) * dphi(1,i);
        dy_dxsi1 += nodes_[connect_(i)] -> getCoordinateValue(1) * dphi(0,i);
        dy_dxsi2 += nodes_[connect_(i)] -> getCoordinateValue(1) * dphi(1,i);        
    };

    //Computing the jacobian determinant
    djac_ = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;
    //djac_ = fabs(djac_);

    //Computing Jacobian inverse
    ainv_(0,0) =  dy_dxsi2 / djac_;
    ainv_(0,1) = -dy_dxsi1 / djac_;
    ainv_(1,0) = -dx_dxsi2 / djac_;
    ainv_(1,1) =  dx_dxsi1 / djac_;

    return;
    
};

//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSpatialDerivatives(ublas::bounded_vector<double,2>& xsi) {
    
    typename QuadShapeFunction<2>::ValueDDeriv ddphi;
    typename QuadShapeFunction<2>::ValueDeriv  dphi;
    
    shapeQuad.evaluateGradient(xsi,dphi);
    shapeQuad.evaluateHessian(xsi,ddphi);
    
    dphi_dx.clear();
    ddphi_dx.clear();

    //Quadratic shape functions spatial first derivatives
    noalias(dphi_dx) = prod(ainv_,dphi);

    // Computing second derivatives
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    double dx_dxsi11 = 0.;
    double dx_dxsi22 = 0.;
    double dx_dxsi12 = 0.;
    double dy_dxsi11 = 0.;
    double dy_dxsi22 = 0.;
    double dy_dxsi12 = 0.;

    for (int i=0; i<6; i++){
        dx_dxsi1 += nodes_[connect_(i)] -> getCoordinateValue(0) * dphi(0,i);
        dx_dxsi2 += nodes_[connect_(i)] -> getCoordinateValue(0) * dphi(1,i);
        dy_dxsi1 += nodes_[connect_(i)] -> getCoordinateValue(1) * dphi(0,i);
        dy_dxsi2 += nodes_[connect_(i)] -> getCoordinateValue(1) * dphi(1,i);

        dx_dxsi11 += nodes_[connect_(i)] -> getCoordinateValue(0) * ddphi(0,0)(i);
        dx_dxsi22 += nodes_[connect_(i)] -> getCoordinateValue(0) * ddphi(1,1)(i);
        dx_dxsi12 += nodes_[connect_(i)] -> getCoordinateValue(0) * ddphi(0,1)(i);
        dy_dxsi11 += nodes_[connect_(i)] -> getCoordinateValue(1) * ddphi(0,0)(i);
        dy_dxsi22 += nodes_[connect_(i)] -> getCoordinateValue(1) * ddphi(1,1)(i);
        dy_dxsi12 += nodes_[connect_(i)] -> getCoordinateValue(1) * ddphi(0,1)(i);        
    };        

    ublas::bounded_matrix<double, 3, 3> invM;
    ublas::bounded_vector<double, 3> vecM, resM;
    
    double a = dx_dxsi1 * dx_dxsi1;
    double b = dy_dxsi1 * dy_dxsi1;
    double c = 2. * dx_dxsi1 * dy_dxsi1;
    double d = dx_dxsi2 * dx_dxsi2;
    double e = dy_dxsi2 * dy_dxsi2;
    double f = 2. * dx_dxsi2 * dy_dxsi2;
    double g = dx_dxsi1 * dx_dxsi2;
    double h = dy_dxsi1 * dy_dxsi2;
    double i = dx_dxsi1 * dy_dxsi2 + dx_dxsi2 * dy_dxsi1;
    
    // matrix:
    // |a b c|
    // |d e f|
    // |g h i|
    
    double det = 1. / (a * e * i - a * f * h - b * d * i + 
                       b * f * g + c * d * h - c * e * g);
    
    invM(0,0) = (e * i - f * h) * det;
    invM(0,1) = (c * h - b * i) * det;
    invM(0,2) = (b * f - c * e) * det;
    invM(1,0) = (f * g - d * i) * det;
    invM(1,1) = (a * i - c * g) * det;
    invM(1,2) = (c * d - a * f) * det;
    invM(2,0) = (d * h - e * g) * det;
    invM(2,1) = (b * g - a * h) * det;
    invM(2,2) = (a * e - b * d) * det;
    
    for (int j = 0; j < 6; j++){
        vecM(0) = ddphi(0,0)(j)
            - dphi_dx(0,j) * dx_dxsi11 - dphi_dx(1,j) * dy_dxsi11;
        vecM(1) = ddphi(1,1)(j)
            - dphi_dx(0,j) * dx_dxsi22 - dphi_dx(1,j) * dy_dxsi22;
        vecM(2) = ddphi(0,1)(j)
            - dphi_dx(0,j) * dx_dxsi12 - dphi_dx(1,j) * dy_dxsi12;

        resM = prod(invM,vecM);
        
        ddphi_dx(0,0)(j) = resM(0);
        ddphi_dx(1,1)(j) = resM(1);
        ddphi_dx(0,1)(j) = resM(2);
        ddphi_dx(1,0)(j) = resM(2);      
    };

    //ddphi_dx.clear();
    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getVelAndDerivatives() {

    u_ = 0.;
    v_ = 0.;

    x_ = 0.;
    y_ = 0.;

    ax_ = 0.;
    ay_ = 0.;

    uPrev_ = 0.;
    vPrev_ = 0.;

    pPrev_ = 0.;

    umesh_ = 0.;
    vmesh_ = 0.;

    p_ = 0.;

    dpPrev_dx = 0.;
    dpPrev_dy = 0.;

    du_dx = 0.;
    du_dy = 0.;
    dv_dx = 0.;
    dv_dy = 0.;

    duprev_dx = 0.;
    duprev_dy = 0.;
    dvprev_dx = 0.;
    dvprev_dy = 0.;

    du_dxx = 0.;
    du_dyy = 0.;
    du_dxy = 0.;
    dv_dxx = 0.;
    dv_dyy = 0.;
    dv_dxy = 0.;

    dp_dx = 0.;
    dp_dy = 0.;

    uInt_ = 0.;
    vInt_ = 0.;
    axInt_ = 0.;
    ayInt_ = 0.;

    dax_dx = 0.;
    dax_dy = 0.;
    day_dx = 0.;
    day_dy = 0.;

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        u_ += nodes_[connect_(i)] -> getVelocity(0) * phi_(i);
        v_ += nodes_[connect_(i)] -> getVelocity(1) * phi_(i);

        uInt_ += nodes_[connect_(i)] -> getIntermediateVelocity(0) * phi_(i);
        vInt_ += nodes_[connect_(i)] -> getIntermediateVelocity(1) * phi_(i);

        axInt_ += nodes_[connect_(i)] -> getIntermediateAcceleration(0) * phi_(i);
        ayInt_ += nodes_[connect_(i)] -> getIntermediateAcceleration(1) * phi_(i);

        x_ += nodes_[connect_(i)] -> getCoordinateValue(0) * phi_(i);
        y_ += nodes_[connect_(i)] -> getCoordinateValue(1) * phi_(i);

        uPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(0) * phi_(i);
        vPrev_ += nodes_[connect_(i)] -> getPreviousVelocity(1) * phi_(i);

        // ax_ += nodes_[connect_(i)] -> getAcceleration(0) * phi_(i);
        // ay_ += nodes_[connect_(i)] -> getAcceleration(1) * phi_(i);
        
        du_dx += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(0,i);
        du_dy += nodes_[connect_(i)] -> getVelocity(0) * dphi_dx(1,i);
        dv_dx += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(0,i);
        dv_dy += nodes_[connect_(i)] -> getVelocity(1) * dphi_dx(1,i);

        duprev_dx += nodes_[connect_(i)] -> getPreviousVelocity(0) * dphi_dx(0,i);
        duprev_dy += nodes_[connect_(i)] -> getPreviousVelocity(0) * dphi_dx(1,i);
        dvprev_dx += nodes_[connect_(i)] -> getPreviousVelocity(1) * dphi_dx(0,i);
        dvprev_dy += nodes_[connect_(i)] -> getPreviousVelocity(1) * dphi_dx(1,i);

        du_dxx += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(0,0)(i);
        du_dyy += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(1,1)(i);
        du_dxy += nodes_[connect_(i)] -> getVelocity(0) * ddphi_dx(0,1)(i);
        dv_dxx += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(0,0)(i);
        dv_dyy += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(1,1)(i);
        dv_dxy += nodes_[connect_(i)] -> getVelocity(1) * ddphi_dx(0,1)(i);

        p_ += nodes_[connect_(i)] -> getPressure() * phi_(i);
        
        pPrev_ += nodes_[connect_(i)] -> getPreviousPressure() * phi_(i);

        dp_dx += nodes_[connect_(i)] -> getPressure() * dphi_dx(0,i);
        dp_dy += nodes_[connect_(i)] -> getPressure() * dphi_dx(1,i);

        dpPrev_dx += nodes_[connect_(i)] -> getPreviousPressure() * dphi_dx(0,i);
        dpPrev_dy += nodes_[connect_(i)] -> getPreviousPressure() * dphi_dx(1,i);

        umesh_ += nodes_[connect_(i)] -> getMeshVelocity(0) * phi_(i);
        vmesh_ += nodes_[connect_(i)] -> getMeshVelocity(1) * phi_(i);
        
        dumesh_dx += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(0,i);
        dumesh_dy += nodes_[connect_(i)] -> getMeshVelocity(0) * dphi_dx(1,i);
        dvmesh_dx += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(0,i);
        dvmesh_dy += nodes_[connect_(i)] -> getMeshVelocity(1) * dphi_dx(1,i);

        dax_dx += (nodes_[connect_(i)]->getVelocity(0) - nodes_[connect_(i)]->getPreviousVelocity(0)) / dTime_ * dphi_dx(0,i);
        dax_dy += (nodes_[connect_(i)]->getVelocity(0) - nodes_[connect_(i)]->getPreviousVelocity(0)) / dTime_ * dphi_dx(1,i);
        day_dx += (nodes_[connect_(i)]->getVelocity(1) - nodes_[connect_(i)]->getPreviousVelocity(1)) / dTime_ * dphi_dx(0,i);
        day_dy += (nodes_[connect_(i)]->getVelocity(1) - nodes_[connect_(i)]->getPreviousVelocity(1)) / dTime_ * dphi_dx(1,i);

    };  

    ax_ = (u_ - uPrev_) / dTime_;
    ay_ = (v_ - vPrev_) / dTime_;

    return;
};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeDragAndLiftForces() {
    
    ublas::bounded_vector<double, 3> nodesb_; 
    if(sideBoundary_ == 0){
        nodesb_(0) = connect_(1); 
        nodesb_(1) = connect_(4); 
        nodesb_(2) = connect_(2); 
        for (int i=0; i<2; i++){
            localNodesBoundary_(0,i) = nodes_[connect_(1)] -> getCoordinateValue(i);
            localNodesBoundary_(1,i) = nodes_[connect_(4)] -> getCoordinateValue(i);
            localNodesBoundary_(2,i) = nodes_[connect_(2)] -> getCoordinateValue(i);
        };
    }else{
        if(sideBoundary_ == 1){
            nodesb_(0) = connect_(2); 
            nodesb_(1) = connect_(5); 
            nodesb_(2) = connect_(0); 
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = nodes_[connect_(2)] -> getCoordinateValue(i);
                localNodesBoundary_(1,i) = nodes_[connect_(5)] -> getCoordinateValue(i);
                localNodesBoundary_(2,i) = nodes_[connect_(0)] -> getCoordinateValue(i);
            };
        }else{
            nodesb_(0) = connect_(0);
            nodesb_(1) = connect_(3);
            nodesb_(2) = connect_(1);
            for (int i=0; i<2; i++){
                localNodesBoundary_(0,i) = nodes_[connect_(0)] -> getCoordinateValue(i);
                localNodesBoundary_(1,i) = nodes_[connect_(3)] -> getCoordinateValue(i);
                localNodesBoundary_(2,i) = nodes_[connect_(1)] -> getCoordinateValue(i);
            };
        };        
    };

    BoundaryQuad           bQuad;     //Boundary Integration Quadrature
    std::pair<BoundaryQuad::PointCoord,BoundaryQuad::PointWeight> gaussQuad;
    DimVector                                  n_vector;
    DimMatrix                                  shearStress;
    ublas::identity_matrix<double>             ident(2);
    ublas::bounded_vector<double,2>            load_friction;
    ublas::bounded_vector<double,2>            load_pressure;
    
    typename QuadShapeFunction<2>::Coords xsi;

    n_vector.clear();
    gaussQuad.first.clear();
    gaussQuad.second.clear();
    shearStress.clear();

    load_friction.clear();
    load_pressure.clear();

    gaussQuad = bQuad.GaussQuadrature();

    double moment = 0.;
    double per = 0.;
    
    int index = 0;
    for(typename BoundaryQuad::QuadratureListIt it = bQuad.begin(); 
        it != bQuad.end(); it++){
        
        double xsiB = gaussQuad.first(index);
        double weightB = gaussQuad.second(index);

        if(sideBoundary_ == 2){
            xsi(0) = (-xsiB + 1.) / 2.;
            xsi(1) = 0.;
        };
        if(sideBoundary_ == 1){
            xsi(1) = (xsiB + 1.) / 2.;
            xsi(0) = 0.;
        };
        if(sideBoundary_ == 0){
            xsi(0) = (xsiB + 1.) / 2.;
            xsi(1) = 1. - xsi(0);
        };

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();        

        phib_ = shapeBound.getShapeFunction(gaussQuad.first(index));
        dphib_ = shapeBound.getShapeFunctionDerivative(gaussQuad.first(index));

        double Tx=0.; double Ty = 0.;
        
        for (int i=0; i<3; i++){
            Tx += localNodesBoundary_(i,0) * dphib_(i);
            Ty += localNodesBoundary_(i,1) * dphib_(i);
        };

        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector(0) =  Ty / jacb_;
        n_vector(1) = -Tx / jacb_;

        shearStress(0,0) = 2. * visc_ * du_dx;
        shearStress(0,1) = visc_ * (du_dy + dv_dx);
        shearStress(1,0) = visc_ * (du_dy + dv_dx);
        shearStress(1,1) = 2. * visc_ * dv_dy;
        
        load_pressure += -p_ * prod(ident,n_vector) * jacb_ * weightB;
        load_friction += prod(shearStress,n_vector) * jacb_ * weightB;

        moment += ((-p_ + shearStress(0,0) + shearStress(1,0)) * y_ 
                 -(-p_ + shearStress(0,1) + shearStress(1,1)) * (x_ - 0.248792267683901))
                 * jacb_ * weightB;
        per += jacb_ * weightB;
        
        index++;
    };

    perimeter = per;
    pitchingMoment = moment;

    pressureDragForce = load_pressure(0);
    pressureLiftForce = load_pressure(1);
    
    frictionDragForce = load_friction(0);
    frictionLiftForce = load_friction(1);

    dragForce = pressureDragForce + frictionDragForce;
    liftForce = pressureLiftForce + frictionLiftForce;

  
    return;
};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getParameterSUPG() {

    DimVector   r, s;

    tSUPG_ = 0.;
    tSUGN1_ = 0.;
    tSUGN2_ = 0.;
    tSUGN3_ = 0.;
    hRGN_ = 0.;
    double hUGN_ = 0.;
    
    r.clear(); s.clear();

    double u__ = 0.;
    double v__ = 0.;

    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);

        double uma = nodes_[connect_(i)] -> getMeshVelocity(0);
        double vma = nodes_[connect_(i)] -> getMeshVelocity(1);
        
        ua -= uma;
        va -= vma;
        
        u__ += ua * phi_(i);
        v__ += va * phi_(i);
    };

    double uNorm = sqrt(u__ * u__ + v__ * v__);
    if(uNorm > 1.e-10){
        s(0) = u__ / uNorm;
        s(1) = v__ / uNorm;
    }else{
        s(0) = 1. / sqrt(2.);
        s(1) = 1. / sqrt(2.);
    };

    for (int i = 0; i < 6; i++){
        double ua = nodes_[connect_(i)] -> getVelocity(0);
        double va = nodes_[connect_(i)] -> getVelocity(1);

        double uma = nodes_[connect_(i)] -> getMeshVelocity(0);
        double vma = nodes_[connect_(i)] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;
        
        r(0) += sqrt(ua * ua + va * va) * dphi_dx(0,i);
        r(1) += sqrt(ua * ua + va * va) * dphi_dx(1,i);
    };

    double rNorm = norm_2(r);

    if (rNorm >= 1.e-10){
        r /= rNorm;
    }else{
        r(0) = 1. / sqrt(2.);
        r(1) = 1. / sqrt(2.);
    };
    
    for (int i = 0; i < 6; i++){
        hRGN_ += fabs(r(0) * dphi_dx(0,i) + r(1) * dphi_dx(1,i));
        hUGN_ += fabs(s(0) * dphi_dx(0,i) + s(1) * dphi_dx(1,i));        
    };

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };    

    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

   
    if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;

 
    //Computing tSUPG parameter
    tSUPG_ = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_));
    //+  1. / (tSUGN3_ * tSUGN3_));
   
    tPSPG_ = tSUPG_;
    tLSIC_ = tSUPG_ * uNorm * uNorm;
 
    return;
};

//------------------------------------------------------------------------------
//-----------------------------ELEMENT LOCAL MATRIX-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemMatrix(int index){
    
    // tSUPG_ = 0.;
    //tPSPG_ = 0.;
    // tLSIC_ = 0.;
 
    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    double duna_dx = timeScheme_ * du_dx + (1. - timeScheme_) * duprev_dx;
    double duna_dy = timeScheme_ * du_dy + (1. - timeScheme_) * duprev_dy;
    double dvna_dx = timeScheme_ * dv_dx + (1. - timeScheme_) * dvprev_dx;
    double dvna_dy = timeScheme_ * dv_dy + (1. - timeScheme_) * dvprev_dy;


    double acelx = (u_ - uPrev_) / dTime_;
    double acely = (v_ - vPrev_) / dTime_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double wSUPGi = (una_ - umesh_) * dphi_dx(0,i) +
                (vna_ - vmesh_) * dphi_dx(1,i);
            double wSUPGj = (una_ - umesh_) * dphi_dx(0,j) +      
                (vna_ - vmesh_) * dphi_dx(1,j);
            
            //Mass matrix (use for both directions)
            double mM =  phi_(i) * phi_(j) * dens_ + 
                wSUPGi * phi_(j) * tSUPG_ * dens_;
            double sMxx = dphi_dx(0,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMxy = dphi_dx(1,i) * phi_(j) * acelx * tSUPG_ * dens_;
            double sMyx = dphi_dx(0,i) * phi_(j) * acely * tSUPG_ * dens_;
            double sMyy = dphi_dx(1,i) * phi_(j) * acely * tSUPG_ * dens_;
                       
            //Difusion matrix (viscosity)
            double Kxx = 2. * dphi_dx(0,i) * dphi_dx(0,j) * visc_ + 
                dphi_dx(1,i) * dphi_dx(1,j) * visc_;
            double Kxy = dphi_dx(1,i) * dphi_dx(0,j) * visc_;
            double Kyx = dphi_dx(0,i) * dphi_dx(1,j) * visc_;
            double Kyy = 2. * dphi_dx(1,i) * dphi_dx(1,j) * visc_ + 
                dphi_dx(0,i) * dphi_dx(0,j) * visc_;
            
            //Convection matrix
            double Cxx = (dphi_dx(0,j) * (una_ - umesh_) + 
                          dphi_dx(1,j) * (vna_ - vmesh_)) * phi_(i) * dens_
                + wSUPGi * wSUPGj * tSUPG_ * dens_;
            
            double Cyy  = Cxx;
            
            double Cuu = phi_(i) * duna_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * duna_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (una_ * duna_dx + vna_ * duna_dy) * dens_;
            
            double Cuv = phi_(i) * duna_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * duna_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (una_ * duna_dx + vna_ * duna_dy) * dens_;
            
            double Cvu = phi_(i) * dvna_dx * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dvna_dx * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(0,i) * phi_(j) * 
                (una_ * dvna_dx + vna_ * dvna_dy) * dens_;
            
            double Cvv = phi_(i) * dvna_dy * phi_(j) * dens_ +
                tSUPG_ * wSUPGi * dvna_dy * phi_(j) * dens_ + 
                tSUPG_ * dphi_dx(1,i) * phi_(j) * 
                (una_ * dvna_dx + vna_ * dvna_dy) * dens_;
            
            double Quu = dphi_dx(0,i) * dp_dx * phi_(j) * tSUPG_;
            double Quv = dphi_dx(1,i) * dp_dx * phi_(j) * tSUPG_;
            double Qvu = dphi_dx(0,i) * dp_dy * phi_(j) * tSUPG_;
            double Qvv = dphi_dx(1,i) * dp_dy * phi_(j) * tSUPG_;
            
            double KLSxx = dphi_dx(0,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSxy = dphi_dx(0,i) * dphi_dx(1,j) * tLSIC_ * dens_;
            double KLSyx = dphi_dx(1,i) * dphi_dx(0,j) * tLSIC_ * dens_;
            double KLSyy = dphi_dx(1,i) * dphi_dx(1,j) * tLSIC_ * dens_;

            double Sxx = + wSUPGi * (2. * ddphi_dx(0,0)(j) + 
                                     ddphi_dx(1,1)(j)) * tSUPG_ * visc_*0;
            double Sxy = + wSUPGi * ddphi_dx(0,1)(j)
                * tSUPG_ * visc_*0;
            double Syx = + wSUPGi * ddphi_dx(0,1)(j)
                * tSUPG_ * visc_*0;
            double Syy = + wSUPGi * (2. * ddphi_dx(1,1)(j) + 
                                     ddphi_dx(0,0)(j)) * tSUPG_ * visc_*0;
            
            
            jacobianNRMatrix(2*i  ,2*j  ) += (mM + timeScheme_ * dTime_ * 
                                              (sMxx + Kxx + KLSxx
                                               + Cxx + Cuu + Quu + Sxx))
                * weight_ * djac_;

            jacobianNRMatrix(2*i+1,2*j+1) += (mM + timeScheme_ * dTime_ *
                                              (sMyy + Kyy + KLSyy
                                               + Cyy + Cvv + Qvv + Syy))
                * weight_ * djac_;

            jacobianNRMatrix(2*i  ,2*j+1) += (timeScheme_ * dTime_ * 
                (Kxy + sMxy + Cuv + Quv + Sxy + KLSxy))
                * weight_ * djac_;

            jacobianNRMatrix(2*i+1,2*j  ) += (timeScheme_ * dTime_ *
                (Kyx + sMyx + Cvu + Qvu + Syx + KLSyx))
                * weight_ * djac_; 

            
            //SINAL DA PARCELA QUE MULTIPLICA O TSUPG ESTA
            //COM SINAL TROCADO NA FORMULAÇAO DO TEZDUYAR
            //multipy pressure direction x
            double QSUPGx = - (dphi_dx(0,i) * phi_(j)) + 
                ((dphi_dx(0,i) * (una_ - umesh_) + 
                  dphi_dx(1,i) * (vna_ - vmesh_)) * 
                 dphi_dx(0,j) * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx(1,i) * phi_(j)) +
                ((dphi_dx(0,i) * (una_ - umesh_) + 
                  dphi_dx(1,i) * (vna_ - vmesh_)) *   
                 dphi_dx(1,j) * tSUPG_);
            //multiply velocity direction x
            double Qx = dphi_dx(0,i) * phi_(j);
            //multiply velocity direction y
            double Qy = dphi_dx(1,i) * phi_(j);
            
            
            jacobianNRMatrix(12+j,2*i  ) += timeScheme_ * Qx * dTime_ * weight_ * djac_;
            jacobianNRMatrix(12+j,2*i+1) += timeScheme_ * Qy * dTime_ * weight_ * djac_;
            jacobianNRMatrix(2*i  ,12+j) += QSUPGx * dTime_ * weight_ * djac_;
            jacobianNRMatrix(2*i+1,12+j) += QSUPGy * dTime_ * weight_ * djac_;

            double Hx = dphi_dx(0,i) * phi_(j) * tPSPG_;
            double Hy = dphi_dx(1,i) * phi_(j) * tPSPG_;
            
            double Gx = dphi_dx(0,i) * wSUPGj * tPSPG_;
            double Gy = dphi_dx(1,i) * wSUPGj * tPSPG_;
            
            double Tx = -(dphi_dx(0,j) * (2. * ddphi_dx(0,0)(i) + 
                                         ddphi_dx(1,1)(i)) + 
                         dphi_dx(1,j) * ddphi_dx(0,1)(i))
                * tPSPG_ * visc_ / dens_*0;
            double Ty = -(dphi_dx(1,j) * (2. * ddphi_dx(1,1)(i) + 
                                         ddphi_dx(0,0)(i)) +
                         dphi_dx(0,j) * ddphi_dx(0,1)(i))
                * tPSPG_ * visc_ / dens_*0;

            double Guu = (dphi_dx(0,i) * duna_dx * phi_(j) + 
                          dphi_dx(1,i) * dvna_dx * phi_(j)) * tPSPG_;
            double Gvv = (dphi_dx(0,i) * duna_dy * phi_(j) + 
                          dphi_dx(1,i) * dvna_dy * phi_(j)) * tPSPG_;
            
            double Q = (dphi_dx(0,i) * dphi_dx(0,j) + 
                        dphi_dx(1,i) * dphi_dx(1,j)) * tPSPG_ / (dens_);
            
            //Hx = 0; Hy = 0; Gx = 0; Gy = 0; Guu = 0; Gvv = 0;


            jacobianNRMatrix(12+j,2*i  ) += (Hx + (Gx + Guu + Tx) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_;
            jacobianNRMatrix(12+j,2*i+1) += (Hy + (Gy + Gvv + Ty) * 
                                              timeScheme_ * dTime_)
                * weight_ * djac_;
            jacobianNRMatrix(12+j,12+i) += Q * dTime_ * weight_ * djac_;
        };
    };

    return;
};


//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditions(){

    for (int i = 0; i < 6; i++){

        if ((nodes_[connect_(i)] -> getConstrains(0) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(0) == 3))  {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            rhsVector(2*i  ) = 0.;
        };

        if ((nodes_[connect_(i)] -> getConstrains(1) == 1) ||
            (nodes_[connect_(i)] -> getConstrains(1) == 3)) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            rhsVector(2*i+1) =  0.;
        };
    };


    // typename Nodes::VecLocD x;

    // for (int i = 0; i < 6; i++){
    //     //if(model){
    //     x = nodes_[connect_(i)] -> getCoordinates();
    //     // double dist = sqrt((x(0)-0.5)*(x(0)-0.5) + (x(1)-0.5)*(x(1)-0.5));
    //     // if(dist < 0.001){
    //     if((x(0) > 0.999) && (x(1) > 0.999)){
    //         // std::cout << "AQUI  " << index_ << std::endl;
            
    //         for (int j = 0; j < 18; j++){
    //             jacobianNRMatrix(12+i,j) = 0.;
    //             jacobianNRMatrix(j,12+i) = 0.;
    //         };
    //         jacobianNRMatrix(12+i,12+i) = 1.;
    //         rhsVector(12+i) =  0.;
    //         //};
    //     };
    // };

    return;
};

//------------------------------------------------------------------------------
//--------------------APPLY THE DIRICHLET BOUNDARY CONDITIONS-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditionsLaplace(){

    for (int i = 0; i < 6; i++){

        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i  ,j) = 0.;
                //jacobianNRMatrix(j,2*i  ) = 0.;
            };
            jacobianNRMatrix(2*i  , 2*i  ) = 1.;
            //rhsVector(2*i  ) = 0.;
        };

        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1) {
            for (int j = 0; j < 18; j++){
                jacobianNRMatrix(2*i+1,j) = 0.;
                //jacobianNRMatrix(j,2*i+1) = 0.;
            };
            jacobianNRMatrix(2*i+1, 2*i+1) = 1.;
            //rhsVector(2*i+1) =  0.;
        };
    };
    
    return;
};


//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getResidualVector(int index){
    
    double una_ = timeScheme_ * u_ + (1. - timeScheme_) * uPrev_;
    double vna_ = timeScheme_ * v_ + (1. - timeScheme_) * vPrev_;

    double duna_dx = timeScheme_ * du_dx + (1. - timeScheme_) * duprev_dx;
    double duna_dy = timeScheme_ * du_dy + (1. - timeScheme_) * duprev_dy;
    double dvna_dx = timeScheme_ * dv_dx + (1. - timeScheme_) * dvprev_dx;
    double dvna_dy = timeScheme_ * dv_dy + (1. - timeScheme_) * dvprev_dy;
    
    for (int i = 0; i < 6; i++){

        // double rMx = dens_ * ((u_ - uPrev_) / dTime_ + 
        //                       (u_ - umesh_) * du_dx)
        //     - dp_dx + visc_ * (ddu_dx + ddu_dy);

        // double rMy = dens_ * ((v_ - vPrev_) / dTime_ + 
        //                       (v_ - vmesh_) * dv_dx)
        //     - dp_dx + visc_ * (ddv_dx + ddv_dy);
        
        double mx = phi_(i) * (u_ - uPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (u_ - uPrev_) * tSUPG_ * dens_;
        double my = phi_(i) * (v_ - vPrev_) * dens_ + 
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) * 
            (v_ - vPrev_) * tSUPG_ * dens_;
        
        double Kx = 2. * dphi_dx(0,i) * duna_dx * visc_ + 
            dphi_dx(1,i) * duna_dy * visc_ + 
            dphi_dx(1,i) * dvna_dx * visc_;
        double Ky= dphi_dx(0,i) * duna_dy * visc_ + 
            2. * dphi_dx(1,i) * dvna_dy * visc_ + 
            dphi_dx(0,i) * dvna_dx * visc_;            

        double KLSx = dphi_dx(0,i) * (duna_dx + dvna_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx(1,i) * (duna_dx + dvna_dy) * tLSIC_ * dens_;

        double Cx = (duna_dx * (una_ - umesh_) + 
                     duna_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * duna_dx + (vna_ - vmesh_) * duna_dy) * tSUPG_ *dens_;
        double Cy = (dvna_dx * (una_ - umesh_) + 
                     dvna_dy * (vna_ - vmesh_)) * phi_(i) * dens_ +
            ((una_ - umesh_) * dphi_dx(0,i) + (vna_ - vmesh_) * dphi_dx(1,i)) *
            ((una_ - umesh_) * dvna_dx + (vna_ - vmesh_) * dvna_dy) * tSUPG_ *dens_;
 
        double Px = - (dphi_dx(0,i) * p_) + ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dx * tSUPG_);
        double Py = - (dphi_dx(1,i) * p_) + ((dphi_dx(0,i) * (una_ - umesh_) +
                                              dphi_dx(1,i) * (vna_ - vmesh_))
                                             * dp_dy * tSUPG_);
           
        double Q = ((duna_dx + dvna_dy) * phi_(i)) +
            (dphi_dx(0,i) * dp_dx + dphi_dx(1,i) * dp_dy) * tPSPG_ / dens_ +
            dphi_dx(0,i) * (u_ - uPrev_) / dTime_ * tPSPG_ +
            dphi_dx(1,i) * (v_ - vPrev_) / dTime_ * tPSPG_ +
            dphi_dx(0,i) * ((una_ - umesh_) * duna_dx +
                            (vna_ - vmesh_) * duna_dy) * tPSPG_ +
            dphi_dx(1,i) * ((una_ - umesh_) * dvna_dx +
                            (vna_ - vmesh_) * dvna_dy) * tPSPG_;

        double T = -(dphi_dx(0,i) * (2. * du_dxx + du_dyy + dv_dxy) +
                    dphi_dx(1,i) * (2. * dv_dyy + dv_dxx + du_dxy))
            * tPSPG_ * visc_ / dens_*0;

        double dAx = 0.;
        double dAy = 0.;

        double Sx = - ((2. * du_dxx + du_dyy + dv_dxy) * 
                       ((una_ - umesh_) * dphi_dx(0,i) + 
                        (vna_ - vmesh_) * dphi_dx(1,i)))
            * tSUPG_ * visc_*0;
        
        double Sy = - ((2. * dv_dyy + dv_dxx + du_dxy) * 
                     ((una_ - umesh_) * dphi_dx(0,i) + 
                      (vna_ - vmesh_) * dphi_dx(1,i)))
            * tSUPG_ * visc_*0;
        

        rhsVector(2*i  ) += (-mx + (-Kx - Px - Cx - dAx - Sx) * dTime_ 
                             - KLSx * dTime_)
            * weight_ * djac_;
        rhsVector(2*i+1) += (-my + (-Ky - Py - Cy - dAy - Sy) * dTime_ 
                             - KLSy * dTime_)
            * weight_ * djac_;
        rhsVector(12+i) += (-Q - T) * dTime_ 
            * weight_ * djac_;
                             
    };

   
    return;
};

//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------

template<>
void Element<2>::getResidualVectorLaplace(){

    rhsVector.clear();
    
    LocalVector U_;
    typename Nodes::VecLocD x_up, x_;

    U_.clear();

    for (int i=0; i<6; i++){

        x_up = nodes_[connect_(i)] -> getUpdatedCoordinates();        
        x_ = nodes_[connect_(i)] -> getCoordinates();

        if (nodes_[connect_(i)] -> getConstrainsLaplace(0) == 1){
            U_(2*i  ) = x_up(0) - x_(0);
        };
        if (nodes_[connect_(i)] -> getConstrainsLaplace(1) == 1){
            U_(2*i+1) = x_up(1) - x_(1);
        };
    };
    
    noalias(rhsVector) += U_;

    return;
};

//------------------------------------------------------------------------------
//---------------------------ELEMENT LAPLACIAN MATRIX---------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getElemLaplMatrix(){

     for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){        
            jacobianNRMatrix(2*i  ,2*j  ) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                              dphi_dx(1,i) * dphi_dx(1,j)) 
                * weight_ * djac_ * meshMovingParameter;
            jacobianNRMatrix(2*i+1,2*j+1) += (dphi_dx(0,i) * dphi_dx(0,j) +
                                              dphi_dx(1,i) * dphi_dx(1,j)) 
                * weight_ * djac_ * meshMovingParameter;
        };
    };
     
    return;
};

//------------------------------------------------------------------------------
//------------------------COMPUTES THE VORTICITY FIELD--------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeVorticity(){
    
    ublas::bounded_matrix <double, 6, 6>      least_squares;
    
    typename QuadShapeFunction<2>::Coords xsi;
    ublas::bounded_vector <double, 6>         nodal_values,results;
    
    int in = 0;
    nodal_values.clear();
    least_squares.clear();
    results.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(in,0);
        xsi(1) = nQuad.PointList(in,1);
                
        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix(xsi);
        
        //Computes spatial derivatives
        getSpatialDerivatives(xsi);
        
        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        for (int i=0; i<6; i++){
            nodal_values(i) += (-du_dy + dv_dx) * phi_(i);
        };
                        
        in++;
    };    

    least_squares.clear();
    
    least_squares(0,0) =  4.825;
    least_squares(0,1) =  1.025;
    least_squares(0,2) =  1.025;
    least_squares(0,3) = -0.275480769230769;  
    least_squares(0,4) =  0.712019230769232;
    least_squares(0,5) = -0.275480769230769;

    least_squares(1,0) =  least_squares(0,1);
    least_squares(1,1) =  4.825;
    least_squares(1,2) =  1.025;
    least_squares(1,3) = -0.275480769230769;  
    least_squares(1,4) = -0.275480769230769;
    least_squares(1,5) =  0.712019230769232;    

    least_squares(2,0) =  least_squares(0,2);
    least_squares(2,1) =  least_squares(1,2);
    least_squares(2,2) =  4.825;
    least_squares(2,3) =  0.712019230769232;  
    least_squares(2,4) = -0.275480769230769;
    least_squares(2,5) = -0.275480769230769;

    least_squares(3,0) =  least_squares(0,3);
    least_squares(3,1) =  least_squares(1,3);
    least_squares(3,2) =  least_squares(2,3);
    least_squares(3,3) =  1.304890902366860;  
    least_squares(3,4) = -0.432609097633136;
    least_squares(3,5) = -0.432609097633136;

    least_squares(4,0) =  least_squares(0,4);
    least_squares(4,1) =  least_squares(1,4);
    least_squares(4,2) =  least_squares(2,4);
    least_squares(4,3) =  least_squares(3,4);  
    least_squares(4,4) =  1.304890902366860;
    least_squares(4,5) = -0.432609097633136;

    least_squares(5,0) =  least_squares(0,5);
    least_squares(5,1) =  least_squares(1,5);
    least_squares(5,2) =  least_squares(2,5);
    least_squares(5,3) =  least_squares(3,5);  
    least_squares(5,4) =  least_squares(4,5);
    least_squares(5,5) =  1.304890902366860;

    results = prod(least_squares,nodal_values);

    for (int i=0; i<6; i++){
        nodes_[connect_(i)] -> incrementVorticity(results(i));
    };
    
    return;
};



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getTransientNavierStokes(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    jacobianNRMatrix.clear();
    rhsVector.clear();

    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){

        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives();

        //Compute Stabilization Parameters
        getParameterSUPG();

        //Computes the element matrix
        getElemMatrix(index);

        //Computes the RHS vector
        getResidualVector(index); 

        index++;        
    };  
    
    // if (index_ == 0){
    //     std::cout << "Element Matrix " << std::endl;
    //     for (int i = 0; i < 18; ++i) {
    //         for (int j = 0; j < 18; ++j){
    //             std::cout << jacobianNRMatrix(i,j) << " " ;         
    //         }
    //         std::cout << std::endl; 
    //     }
    // };


    //Apply boundary conditions
    setBoundaryConditions();

    return;
};

//------------------------------------------------------------------------------
//----------------------------STEADY LAPLACE PROBEM-----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSteadyLaplace(){

    typename QuadShapeFunction<2>::Coords xsi;
    int index = 0;

    jacobianNRMatrix.clear();
    
    for(typename NormalQuad::QuadratureListIt it = nQuad.begin(); 
        it != nQuad.end(); it++){
        
        //Defines the integration points adimentional coordinates
        xsi(0) = nQuad.PointList(index,0);
        xsi(1) = nQuad.PointList(index,1);

        //Computes the velocity shape functions
        shapeQuad.evaluate(xsi,phi_);

        //Returns the quadrature integration weight
        weight_ = nQuad.WeightList(index);

        //Computes the jacobian matrix
        getJacobianMatrix(xsi);

        //Computes spatial derivatives
        getSpatialDerivatives(xsi);

        //Computes the element matrix
        getElemLaplMatrix();

        index++;        
    };  
    
    //Computes the RHS vector
    getResidualVectorLaplace();

    //Apply boundary conditions
    setBoundaryConditionsLaplace();

    return;
};



#endif

