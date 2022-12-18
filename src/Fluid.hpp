//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------FLUID-------------------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_H
#define FLUID_H

#include "Element.hpp"
#include "Boundary.hpp"
#include "fluidDomain.h"
#include "StructuralDomain.hpp"

// PETSc libraries
#include <metis.h>
#include <petscksp.h> 

#include <boost/timer.hpp> 
#include <boost/thread.hpp>

/// Mounts the incompressible flow problem
template<int DIM>
class Fluid{
public:
    /// Defines the class Element locally
    typedef Element<DIM> Elements;

    /// Defines the class Node locally
    typedef typename Elements::Nodes  Node;

    /// Defines the class Boundary locally
    typedef Boundary<DIM> Boundaries;

    /// Defines the vector of fluid nodes
    std::vector<Node *>       nodes_;

    /// Defines the vector of fluid elements
    std::vector<Elements *>   elements_;
 
    /// Defines the vector of fluid boundaries mesh nodes
    std::vector<Boundaries *> boundary_;

private:
    //FLUID VARIABLES
    std::string inputFile; //Fluid input file
    int numElem;           //Number of elements in fluid mesh for each thread 
    int numTotalElem;           //Total number of elements in fluid mesh 
    int numNodes;          //Number of nodes in velocity/quadratic mesh
    int numBoundaries;     //Number of fluid boundaries
    int numBoundElems;     //Number of elements in fluid boundaries
    double pressInf;       //Undisturbed pressure 
    double rhoInf;         //Density
    double viscInf;        //Viscosity
    double velocityInf[3]; //Undisturbed velocity
    double fieldForces[3]; //Field forces (constant)
    idx_t* part_elem;      //Fluid Domain Decomposition - Elements
    idx_t* part_nodes;     //Fluid Domain Decomposition - Nodes
    int numTimeSteps;      //Number of Time Steps
    int printFreq;         //Printing frequence of output files
    double dTime;          //Time Step
    double integScheme;    //Time Integration Scheme
    int rank;
    int size;
    bool computeDragAndLift;
    int numberOfLines;
    std::vector<int> dragAndLiftBoundary;
    int iTimeStep;
    Geometry* geometry_;
    double pi = M_PI;

    // FSI variables
    StructuralDomain structure_;
    int initialStructuralTimeStep_ = 0;
    double firstTimeStep_;

    
public:
    int weightFunctionBehavior;
    bool printVelocity;
    bool printPressure;
    bool printVorticity;
    bool printMeshVelocity;
    bool printMeshDisplacement;
    bool printJacobian;
    bool printProcess;

public:
    /// Reads the input file and perform the preprocessing operations
    /// @param Geometry* mesh geometry @param std::string input file 
    /// @param std::string input .msh file @param std::string mirror file
    /// @param bool delete mesh files
    void dataReading(Geometry* geometry, const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles, const std::string& structuralInputFile);

    void readInitialValues(const std::string& inputVel,const std::string& inputPres);


    /// Performs the domain decomposition for parallel processing
    void domainDecompositionMETIS(std::vector<Elements *> &elem_); 

    /// Export the domain decomposition 
    /// @return pair with the elements and nodes domain decompositions
    std::pair<idx_t*,idx_t*> getDomainDecomposition(){
        return std::make_pair(part_elem,part_nodes);};
    
    /// Gets the flag for computing the Drag and Lift coefficients in a 
    /// specific boundary
    /// @return bool flag for computing Drag and Lift coefficients
    bool getComputeDragAndLift() {return computeDragAndLift;};

    /// Gets the number of the boundary for computing the Drag and Lift 
    /// coefficients
    /// @return int boundary number
    int getDragAndLiftBoundary(int index) {return dragAndLiftBoundary[index];};
 
    /// Mounts and solve the transient incompressible flow problem    
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    int solveTransientProblem(int iterNumber,double tolerance);

    /// Mounts and solve the transient incompressible flow problem for moving
    /// domain problems
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    int solveTransientProblemMoving(int iterNumber,double tolerance);

    /// Mounts and solve the steady Laplace problem
    /// @param int maximum number of Newton-Raphson's iterations
    /// @param double tolerance of the Newton-Raphson's process
    int solveSteadyLaplaceProblem(int iterNumber, double tolerance);

    /// Print the results for Paraview post-processing
    /// @param int time step
    void printResults(int step);

    /// Compute and print drag and lift coefficients
    void dragAndLiftCoefficients(std::ofstream& dragLift);

    /// Computes forces on fluid-structure interface
    void computeFSIForces(std::vector<double> &F);

    /// Prints the FSI boundary
    void printFSIBoundary(const int& timestep, std::string output_file);

    /// Prints the structure's center displacements
    void printStructureCenterDisplacements(const int &timeStep, std::ofstream &out);
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::domainDecompositionMETIS(std::vector<Elements *> &elem_) {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);



    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t three = 3;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[(4*dd-2)*numEl];

    MPI_Bcast(&numEl,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(&numNd,1,MPI_INT,0,PETSC_COMM_WORLD);


    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];


    if (rank == 0){
        for (idx_t i = 0; i < numEl+1; i++){
            elem_start[i]=(4*dd-2)*i;
        };
        for (idx_t jel = 0; jel < numEl; jel++){
            typename Elements::Connectivity connec;
            connec=elem_[jel]->getConnectivity();        
            
            for (idx_t i=0; i<(4*dd-2); i++){
            elem_connec[(4*dd-2)*jel+i] = connec(i);
            };
        };

        //Performs the domain decomposition
        METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                                  NULL, NULL, &one, &ssize, NULL, NULL,
                                  &objval, part_elem, part_nodes);

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
        for(int i = 0; i < numElem; i++){
            mirrorData << "process = " << part_elem[i]
                       << ", element = " << i << std::endl;
        };

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
        for(int i = 0; i < numNodes; i++){
            mirrorData << "process = " << part_nodes[i]
                       << ", node = " << i << std::endl;
        };
        

        for (int i = 0; i < size; ++i){
            std::string result;
            std::ostringstream convert;

            convert << i+000;
            result = convert.str();
            std::string s = "mesh"+result+".dat";

            std::fstream mesh(s.c_str(), std::ios_base::out);

            int locElem = std::count(part_elem, part_elem+numElem, i);

            mesh << locElem << std::endl;

            for (int jel = 0; jel < numElem; ++jel){
                if (part_elem[jel] == i){
                    typename Elements::Connectivity connec;
                    connec = elem_[jel]->getConnectivity();
                    mesh << jel << " " << connec(0) << " " << connec(1) << " " << connec(2) << " "
                         << connec(3) << " " << connec(4) << " " << connec(5) << std::endl;
                }
            }

        }
    }

    MPI_Bcast(part_elem,numEl,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(part_nodes,numNd,MPI_INT,0,PETSC_COMM_WORLD);

    return;

};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::printResults(int step) {

    std::string s;
    
    if (step % printFreq == 0){

        std::string result;
        std::ostringstream convert;

        convert << step+100000;
        result = convert.str();
        if (rank == 0) s = "saidaVel"+result+".vtu";

        std::fstream output_v(s.c_str(), std::ios_base::out);

        if (rank == 0){
            output_v << "<?xml version=\"1.0\"?>" << std::endl
                     << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                     << "  <UnstructuredGrid>" << std::endl
                     << "  <Piece NumberOfPoints=\"" << numNodes
                     << "\"  NumberOfCells=\"" << numTotalElem
                     << "\">" << std::endl;

            //WRITE NODAL COORDINATES
            output_v << "    <Points>" << std::endl
                     << "      <DataArray type=\"Float64\" "
                     << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

            for (int i=0; i<numNodes; i++){
                typename Node::VecLocD x;
                x=nodes_[i]->getCoordinates();
                output_v << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
            };
            output_v << "      </DataArray>" << std::endl
                     << "    </Points>" << std::endl;
            
            //WRITE ELEMENT CONNECTIVITY
            output_v << "    <Cells>" << std::endl
                     << "      <DataArray type=\"Int32\" "
                     << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
        }

        int k = 0;
        for (int iElem = 0; iElem < numTotalElem; ++iElem){
            typename Elements::Connectivity connec;

            if (part_elem[iElem] == rank){
                connec = elements_[k]->getConnectivity();
                k++;
            }
            MPI_Bcast(&connec,8,MPI_INT,part_elem[iElem],PETSC_COMM_WORLD);
            if (rank == 0) output_v << connec(0) << " " << connec(1) << " " << connec(2) << " "
                                    << connec(3) << " " << connec(4) << " " << connec(5) << std::endl;

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        if (rank == 0) {
            output_v << "      </DataArray>" << std::endl;
      
            //WRITE OFFSETS IN DATA ARRAY
            output_v << "      <DataArray type=\"Int32\""
                    << " Name=\"offsets\" format=\"ascii\">" << std::endl;
        
            int aux = 0;
            for (int i=0; i<numTotalElem; i++){
                output_v << aux + 6 << std::endl;
                aux += 6;
            };
            output_v << "      </DataArray>" << std::endl;
      
            //WRITE ELEMENT TYPES
            output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                     << "format=\"ascii\">" << std::endl;
        
            for (int i=0; i<numTotalElem; i++){
                output_v << 22 << std::endl;
            };

            output_v << "      </DataArray>" << std::endl
                     << "    </Cells>" << std::endl;

            //WRITE NODAL RESULTS
            output_v << "    <PointData>" << std::endl;

            if (printVelocity){
                output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                        << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                for (int i=0; i<numNodes; i++){
                    output_v << nodes_[i] -> getVelocity(0) << " "             
                             << nodes_[i] -> getVelocity(1) << " " << 0. << std::endl;
                }; 
                output_v << "      </DataArray> " << std::endl;
            };

            if (printMeshVelocity){
                output_v<< "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                        << "Name=\"Mesh Velocity\" format=\"ascii\">" << std::endl;
            
                for (int i=0; i<numNodes; i++){
                    output_v << nodes_[i] -> getMeshVelocity(0) << " "    
                             << nodes_[i] -> getMeshVelocity(1) << " " 
                             << 0. << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (printVorticity){
                output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Vorticity\" format=\"ascii\">" << std::endl;
                for (int i=0; i<numNodes; i++){
                    output_v << nodes_[i] -> getVorticity() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (printPressure){
                output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                for (int i=0; i<numNodes; i++){
                    output_v << 0. << " " << 0. << " " 
                             << nodes_[i] -> getPressure() << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            if (printMeshDisplacement){
                output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                         << "Name=\"Mesh Displacement\" format=\"ascii\">" << std::endl;
                for (int i=0; i<numNodes; i++){       
                    typename Node::VecLocD x, xp;
                    x=nodes_[i]->getCoordinates();
                    xp=nodes_[i]->getInitialCoordinates();
                
                    output_v << x(0)-xp(0) << " " << x(1)-xp(1) << " " 
                            << 0. << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };

            output_v << "    </PointData>" << std::endl; 

            //WRITE ELEMENT RESULTS
            output_v << "    <CellData>" << std::endl;
        };


        if (printProcess){
            if(rank == 0){
                output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Process\" format=\"ascii\">" << std::endl;
            
                for (int i=0; i<numTotalElem; i++){
                    output_v << part_elem[i] << std::endl;
                };
                output_v << "      </DataArray> " << std::endl;
            };
        };

        if (printJacobian){
            if (rank == 0){
                output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                         << "Name=\"Jacobian\" format=\"ascii\">" << std::endl;
            };
            

            int k = 0;
            for (int iElem = 0; iElem < numTotalElem; ++iElem){
                double jac = 0.;

                if (part_elem[iElem] == rank){
                    jac = elements_[k]->getJacobian();
                    k++;
                }

                MPI_Bcast(&jac,1,MPI_DOUBLE,part_elem[iElem],PETSC_COMM_WORLD);
                if (rank == 0) output_v << jac << std::endl;

                MPI_Barrier(PETSC_COMM_WORLD);
            }                

            if(rank == 0) output_v << "      </DataArray> " << std::endl;
        };
    

        if(rank == 0){
            output_v << "    </CellData>" << std::endl; 

            //FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl;
            output_v << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;
         };

    };

return;

};
//------------------------------------------------------------------------------
//----------------------COMPUTES DRAG AND LIFT COEFFICIENTS---------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::dragAndLiftCoefficients(std::ofstream& dragLift){
    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;
    
    for (int jel = 0; jel < numBoundElems; jel++){   
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        
        
        for (int i=0; i<numberOfLines; i++){
            if (boundary_[jel] -> getBoundaryGroup() == dragAndLiftBoundary[i]){
                //std::cout << "AQUI " << numberOfLines<< " " << i << " " << dragAndLiftBoundary[i] << std::endl;
                int iel = boundary_[jel] -> getElement();
                
                for (int j = 0; j < numElem; ++j){
                    if (elements_[j] -> getIndex() == iel){
                        elements_[j] -> computeDragAndLiftForces(structure_.getCenter().getCurrentPositions());
                    
                        pDForce = elements_[j] -> getPressureDragForce();
                        pLForce = elements_[j] -> getPressureLiftForce();
                        fDForce = elements_[j] -> getFrictionDragForce();
                        fLForce = elements_[j] -> getFrictionLiftForce();
                        dForce = elements_[j] -> getDragForce();
                        lForce = elements_[j] -> getLiftForce();
                        pMom += elements_[j] -> getPitchingMoment();
                        per += elements_[j] -> getPerimeter();
                    }   
                }
            };
        };
        
        pressureDragCoefficient += pDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        pressureLiftCoefficient += pLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf[0] * velocityInf[0]);
        
    };



    MPI_Barrier(PETSC_COMM_WORLD);

    double totalDragCoefficient = 0.;
    double totalLiftCoefficient = 0.;
    double totalPressureDragCoefficient = 0.;
    double totalPressureLiftCoefficient = 0.;
    double totalFrictionDragCoefficient = 0.;
    double totalFrictionLiftCoefficient = 0.;
    double totalPitchingMomentCoefficient = 0.;;
    double perimeter = 0.;

    MPI_Allreduce(&dragCoefficient,&totalDragCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&liftCoefficient,&totalLiftCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&pressureDragCoefficient,&totalPressureDragCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&pressureLiftCoefficient,&totalPressureLiftCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&frictionDragCoefficient,&totalFrictionDragCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&frictionLiftCoefficient,&totalFrictionLiftCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&pMom,&pitchingMomentCoefficient,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&per,&perimeter,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    totalPitchingMomentCoefficient = pitchingMomentCoefficient / (rhoInf * velocityInf[0] * velocityInf[0] * perimeter);


    if (rank == 0) {
        const int timeWidth = 15;
        const int numWidth = 15;
        dragLift << std::setprecision(5) << std::scientific;
        if (iTimeStep >= initialStructuralTimeStep_){
            dragLift << std::left << std::setw(timeWidth) << initialStructuralTimeStep_ * firstTimeStep_ + (iTimeStep - initialStructuralTimeStep_) * dTime;
        } else {
            dragLift << std::left << std::setw(timeWidth) << iTimeStep * firstTimeStep_;
        }
        dragLift << std::setw(numWidth) << totalPressureDragCoefficient;
        dragLift << std::setw(numWidth) << totalPressureLiftCoefficient;
        dragLift << std::setw(numWidth) << totalFrictionDragCoefficient;
        dragLift << std::setw(numWidth) << totalFrictionLiftCoefficient;
        dragLift << std::setw(numWidth) << totalDragCoefficient;
        dragLift << std::setw(numWidth) << totalLiftCoefficient;
        dragLift << std::setw(numWidth) << totalPitchingMomentCoefficient;
        dragLift << std::endl;
    }
}
//------------------------------------------------------------------------------
//------------------------COMPUTES FORCES ON FSI INTERFACE----------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::computeFSIForces(std::vector<double> &F){
    double dForce = 0.;
    double lForce = 0.;
    double pMom = 0.;
    
    for (int jel = 0; jel < numBoundElems; jel++){   
        for (int i=0; i<numberOfLines; i++){
            if (boundary_[jel] -> getBoundaryGroup() == dragAndLiftBoundary[i]){
                //std::cout << "AQUI " << numberOfLines<< " " << i << " " << dragAndLiftBoundary[i] << std::endl;
                int iel = boundary_[jel] -> getElement();
                
                for (int j = 0; j < numElem; ++j){
                    if (elements_[j] -> getIndex() == iel){
                        elements_[j] -> computeDragAndLiftForces(structure_.getCenter().getCurrentPositions());
                    
                        dForce += elements_[j] -> getDragForce();
                        lForce += elements_[j] -> getLiftForce();
                        pMom += elements_[j] -> getPitchingMoment();
                    }   
                }
            };
        };
    };

    MPI_Barrier(PETSC_COMM_WORLD);

    double totalDragForce = 0.;
    double totalLiftForce = 0.;
    double totalPitchingMoment = 0.;;

    MPI_Allreduce(&dForce,&totalDragForce,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lForce,&totalLiftForce,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&pMom,&totalPitchingMoment,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    F[0] = totalDragForce;
    F[1] = totalLiftForce;
    F[2] = totalPitchingMoment;
}

template<>
void Fluid<2>::printFSIBoundary(const int& timestep, string output_file){
    
    if (timestep % printFreq == 0){
        if (rank == 0){
            std::stringstream text;
            text << "Results/Uncoupled/Euler/" << output_file << timestep << ".vtu";
            std::ofstream file(text.str());
            file.precision(16);

            //header
            file << "<?xml version=\"1.0\"?>" << "\n"
                << "<VTKFile type=\"UnstructuredGrid\">" << "\n"
                << "  <UnstructuredGrid>" << "\n"
                << "  <Piece NumberOfPoints=\"" << structure_.getBoundary().size()
                << "\"  NumberOfCells=\"" << structure_.getBoundary().size()
                << "\">" << "\n";
            //nodal coordinates
            file << "    <Points>" << "\n"
                << "      <DataArray type=\"Float64\" "
                << "NumberOfComponents=\"3\" format=\"ascii\">" << "\n";
            for (StructuralNode n : structure_.getBoundary()){
                file << n.getCurrentPosition(0) << " " << n.getCurrentPosition(1) << " " << 0.0 << "\n";
            }
            file << "      </DataArray>" << "\n"
                << "    </Points>" << "\n";
            //element connectivity
            file << "    <Cells>" << "\n"
                << "      <DataArray type=\"Int32\" "
                << "Name=\"connectivity\" format=\"ascii\">" << "\n";
            for (int i = 1; i <= structure_.getBoundary().size(); i++){
                file << structure_.getBoundary()[i%structure_.getBoundary().size()].getIndex() << " " << structure_.getBoundary()[i-1].getIndex() << " ";
                file << "\n";
            }

            file << "      </DataArray>" << "\n";
            //offsets
            file << "      <DataArray type=\"Int32\""
                << " Name=\"offsets\" format=\"ascii\">" << "\n";
            int aux = 0;
            for (StructuralNode n : structure_.getBoundary()){
                aux += 2;
                file << aux << "\n";
            }
            file << "      </DataArray>" << "\n";
            //elements type
            file << "      <DataArray type=\"UInt8\" Name=\"types\" "
                << "format=\"ascii\">" << "\n";

            for (StructuralNode n : structure_.getBoundary()){
                file << 3 << "\n";
            }
            file << "      </DataArray>" << "\n"
                << "    </Cells>" << "\n";
            //nodal results
            file << "    <PointData>" <<"\n";
            file << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                << "Name=\"Position\" format=\"ascii\">" << "\n";

            for (StructuralNode n: structure_.getBoundary()){
                std::vector<double> u_ = n.getCurrentPositions();
                file << u_[0] << " "
                    << u_[1]  << " " << 0.0 << "\n";
            }
            file << "      </DataArray> " << "\n";
            file << "    </PointData>" << "\n";
            //elemental results
            file << "    <CellData>" << "\n";

            file << "    </CellData>" << "\n";
            //footnote
            file << "  </Piece>" << "\n"
                << "  </UnstructuredGrid>" << "\n"
                << "</VTKFile>" << "\n";
            file.close();
        }
    }
}

template<>
void Fluid<2>::printStructureCenterDisplacements(const int &timeStep, std::ofstream &out){
    if (rank == 0){
        const int timeWidth = 15;
        const int numWidth = 15;

        out << std::setprecision(8);
        if (timeStep <= initialStructuralTimeStep_){
            out << std::left << std::setw(timeWidth) << timeStep*firstTimeStep_ << " "; 
        } else {
            out << std::left << std::setw(timeWidth) << initialStructuralTimeStep_*firstTimeStep_ + (timeStep-initialStructuralTimeStep_)*dTime;
        }

        out << std::setw(numWidth) << structure_.getCenter().getCurrentPosition(0) - structure_.getCenter().getInitialPosition(0);
            << std::setw(numWidth) << structure_.getCenter().getCurrentPosition(1) - structure_.getCenter().getInitialPosition(1);
            << std::setw(numWidth) << structure_.getCenter().getCurrentPosition(2) - structure_.getCenter().getInitialPosition(2) << endl;

    }
}

//------------------------------------------------------------------------------
//----------------------------READS FLUID INPUT FILE----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::dataReading(Geometry* geometry, const std::string& inputFile, 
                           const std::string& inputMesh, const std::string& mirror,
                           const bool& deleteFiles, const std::string& structuralInputFile){

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if (rank == 0)std::cout << "Reading fluid data from \"" 
                            << inputFile << "\"" << std::endl;

    std::ifstream file(inputMesh);
    std::string line;
    std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  
    
    geometry_ = geometry;


    //defyning the maps that are used to store the elements information
    std::unordered_map<int, std::string> gmshElement = { {1, "line"}, {2, "triangle"}, {3, "quadrilateral"}, {8, "line3"}, {9, "triangle6"}, {10, "quadrilateral9"}, {15, "vertex"}, {16, "quadrilateral8"}, {20, "triangle9"}, {21, "triangle10"}, {26, "line4"}, {36, "quadrilateral16"}, {39, "quadrilateral12"} };
    std::unordered_map<std::string, int> numNodes2 = { {"vertex", 1}, {"line", 2}, {"triangle", 3}, {"quadrilateral", 4}, {"line3", 3}, {"triangle6", 6}, {"quadrilateral8", 8}, {"quadrilateral9", 9}, {"line4", 4}, {"triangle", 9}, {"triangle10", 10}, {"quadrilateral12", 12}, {"quadrilateral16", 16}};
    std::unordered_map<std::string, std::string> supportedElements = { {"triangle", "T3"}, {"triangle6", "T6"}, {"triangle10", "T10"}, {"quadrilateral", "Q4"}, {"quadrilateral8", "Q8"}, {"quadrilateral9", "Q9"}, {"quadrilateral12", "Q12"}, {"quadrilateral16", "Q16"} };
    std::unordered_map<Line*, std::vector< std::vector<int> >> lineElements;


    //Defines input and output files    
    std::ifstream inputData(inputFile.c_str());
    std::ofstream mirrorData(mirror.c_str());

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numTimeSteps >> printFreq;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density and viscosity
    inputData >> rhoInf >> viscInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime >> integScheme;
    firstTimeStep_ = dTime;

    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Drag and lift
    inputData >> computeDragAndLift >> numberOfLines; 
    dragAndLiftBoundary.reserve(numberOfLines);
    for (int i = 0; i < numberOfLines; ++i)
    {
        int aux;
        inputData >> aux; 
        dragAndLiftBoundary.push_back(aux);
    }
    

    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
    for (int i = 0; i < numberOfLines; ++i)
    {
        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
    }



    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Printing results
    inputData >> printVelocity;              getline(inputData,line);
    inputData >> printPressure;              getline(inputData,line);
    inputData >> printVorticity;             getline(inputData,line);
    inputData >> printMeshVelocity;          getline(inputData,line);
    inputData >> printMeshDisplacement;      getline(inputData,line);
    inputData >> printJacobian;              getline(inputData,line);
    inputData >> printProcess;              

    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
    mirrorData << "PrintPressure              = " << printPressure << std::endl;
    mirrorData << "PrintVorticity             = " << printVorticity 
               << std::endl;
    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity
               << std::endl;
    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement
               << std::endl;
    mirrorData << "PrintJacobian              = " << printJacobian << std::endl;
    mirrorData << "PrintProcess               = " << printProcess << std::endl 
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);
 
    int dimension=2;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++READIN MESH+++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++PHYSICAL ENTITIES++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_physical_entities;
    file >> number_physical_entities;
    std::getline(file, line);
    std::unordered_map<int, std::string> physicalEntities;
    physicalEntities.reserve(number_physical_entities);

    for (int i = 0; i < number_physical_entities; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }
    std::getline(file, line); std::getline(file, line);

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++NODES++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    file >> numNodes;
    nodes_.reserve(numNodes);
    std::getline(file, line);
    int index = 0;
    if (rank == 0) std::cout << "Number of Nodes " << " " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++)
    {
        typename Node::VecLocD x;
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        bounded_vector<double,2> coord;
        std::istringstream(tokens[1]) >> x(0);
        std::istringstream(tokens[2]) >> x(1);
        //addNode(i, coord);
         Node *node = new Node(x, index++);
         nodes_.push_back(node);
    }
    std::getline(file, line); std::getline(file, line);

    mirrorData << "Nodal Coordinates " << numNodes << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        typename Node::VecLocD x;
        x = nodes_[i]->getCoordinates();       
        for (int j=0; j<2; j++){
            mirrorData << x(j) << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
        double u[2];
        u[0] = 0.; u[1] = 0.;
        nodes_[i] -> setMeshVelocity(u);
    };


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int number_elements;
    file >> number_elements;
    //elements_.reserve(number_elements);
    std::vector<Elements *>   elementsAux_;
    elementsAux_.reserve(number_elements);

    boundary_.reserve(number_elements/10);
    index = 0;
    std::getline(file, line);
    int cont = 0;

    numBoundElems = 0;
    numElem = 0;

    std::vector<BoundaryCondition*> dirichlet, neumann, glue, FSinterface;
    dirichlet = geometry_->getBoundaryCondition("DIRICHLET"); 
    neumann = geometry_->getBoundaryCondition("NEUMANN"); 
    FSinterface = geometry_->getBoundaryCondition("MOVING");
    glue = geometry_->getBoundaryCondition("GEOMETRY");

    std::set<int> FSIBoundaryNodesNumbers; // Set to identify FSI boundary
    for (int i = 0; i < number_elements; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        std::vector<int> values(tokens.size(), 0);
        for (size_t j = 0; j < tokens.size(); j++)
            std::istringstream(tokens[j]) >> values[j];
        std::string elementType = gmshElement[values[1]];
        int number_nodes_per_element = numNodes2[elementType];
        std::vector<int> elementNodes;
        elementNodes.reserve(number_nodes_per_element);

        for (size_t j = 5 ; j < values.size(); j++)
            elementNodes.push_back(values[j]-1);
 
        std::string name = physicalEntities[values[3]];
        //Adding 2D elements to surfaces
        if (name[0] == 's'){

            if(rank == 0){
                if (supportedElements.find(elementType) == supportedElements.end()){
                    std::cout << elementType << " is not supported.\n";
                    exit(EXIT_FAILURE);
                }

                PlaneSurface* object = geometry_ -> getPlaneSurface(name);
                //int materialIndex = object -> getMaterial() -> getIndex();
                //double thickness = object -> getThickness();
                numElem++;

                typename Elements::Connectivity connect;
                connect.clear();
                for (int j = 0 ; j < 6; j++) connect(j) = elementNodes[j];

                Elements *el = new Elements(index++,connect,nodes_);
                elementsAux_.push_back(el);

                for (int k = 0; k<6; k++){
                    nodes_[connect(k)] -> pushInverseIncidence(index);
                };
            }
        }
        else if (name[0] == 'l')
        {
            Boundaries::BoundConnect connectB;

            connectB(0) = elementNodes[0];
            connectB(1) = elementNodes[1];
            connectB(2) = elementNodes[2];

            int ibound;

            std::string::size_type sz;   // alias of size_t
            ibound = std::stoi (&name[1],nullptr,10);

            int constrain[3];
            double value[3];

            for (int i = 0; i < dirichlet.size(); i++){
                if (name == dirichlet[i] -> getLineName()){
                     if ((dirichlet[i] -> getComponentX()).size() == 0){
                         constrain[0] = 0; value[0] = 0;
                     }else{
                         std::vector<double> c = dirichlet[i] -> getComponentX();
                         constrain[0] = 1;
                         value[0] = c[0];
                     }
                     if ((dirichlet[i] -> getComponentY()).size() == 0){
                         constrain[1] = 0; value[1] = 0;
                     }else{
                         std::vector<double> c = dirichlet[i] -> getComponentY();
                         constrain[1] = 1;
                         value[1] = c[0];
                     }
                }
            }

            for (int i = 0; i < neumann.size(); i++){
                if (name == neumann[i] -> getLineName()){
                    if ((neumann[i] -> getComponentX()).size() == 0){
                        constrain[0] = 0; value[0] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX();
                        constrain[0] = 0;
                        value[0] = c[0];
                    }
                    if ((neumann[i] -> getComponentY()).size() == 0){
                        constrain[1] = 0; value[1] = 0;
                    }else{
                        std::vector<double> c = neumann[i] -> getComponentX(); // AQUI NÃO DEVIA SER Y?
                        constrain[1] = 0;
                        value[1] = c[0];
                    }
                }
            }  
            
            for (int i = 0; i < glue.size(); i++){
                if (name == glue[i] -> getLineName()){
                    if ((glue[i] -> getComponentX()).size() == 0){
                        constrain[0] = 2; value[0] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentX();
                        constrain[0] = 2;
                        value[0] = c[0];
                    }
                    if ((glue[i] -> getComponentY()).size() == 0){
                        constrain[1] = 2; value[1] = 0;
                    }else{
                        std::vector<double> c = glue[i] -> getComponentY();
                        constrain[1] = 2;
                        value[1] = c[0];
                    }//std::cout <<"aqui " << std::endl;
                }
            }              
            for (int i = 0; i < FSinterface.size(); i++){
                if (name == FSinterface[i] -> getLineName()){
                    if ((FSinterface[i] -> getComponentX()).size() == 0){
                        constrain[0] = 3; value[0] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentX();
                        constrain[0] = 3;
                        value[0] = c[0];
                    }
                    if ((FSinterface[i] -> getComponentY()).size() == 0){
                        constrain[1] = 3; value[1] = 0;
                    }else{
                        std::vector<double> c = FSinterface[i] -> getComponentY();
                        constrain[1] = 3;
                        value[1] = c[0];
                    }
                    FSIBoundaryNodesNumbers.insert(connectB(0));
                    FSIBoundaryNodesNumbers.insert(connectB(1));
                    FSIBoundaryNodesNumbers.insert(connectB(2));
                }
            }        
            Boundaries * bound = new Boundaries(connectB, numBoundElems++, constrain, value, ibound);
            // std::cout << "asdasd " << rank << " " << ibound << std::endl;
            boundary_.push_back(bound);           
        }   
    }

    domainDecompositionMETIS(elementsAux_);

    if (rank == 0){
        for (int i = 0; i < numElem; ++i) delete elementsAux_[i];
        elementsAux_.clear();
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    std::string result;
    std::ostringstream convert;

    convert << rank+000;
    result = convert.str();
    std::string s = "mesh"+result+".dat";

    std::ifstream mesh(s.c_str(), std::ios_base::out);

    mesh >> numElem;

    elements_.reserve(numElem);

    //reading element connectivity
    for (int i = 0; i < numElem; i++){
        typename Elements::Connectivity connect;
        connect.clear();
        int ind_ = 0;

        mesh >> ind_ >> connect(0) >> connect(1) >> connect(2) >> connect(3) >> connect(4) >> connect(5);

        Elements *el = new Elements(ind_,connect,nodes_);
        elements_.push_back(el);
    };

    MPI_Barrier(PETSC_COMM_WORLD);

    MPI_Allreduce(&numElem,&numTotalElem,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

    if (rank == 0) std::cout << "Number of elements " << number_elements << " " 
                             << numTotalElem << " " << numBoundElems << std::endl;
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();       
        for (int i=0; i<4*dimension-2; i++){
            mirrorData << connec(i) << " ";
        };
        mirrorData << std::endl;
    };

    //Sets boundary constrains
    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        Boundaries::BoundConnect connectB;
        connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB(0);
        int no2 = connectB(1);
        int no3 = connectB(2);
        if ((boundary_[ibound] -> getConstrain(0) != 2)){
            nodes_[no1] -> setConstrainsLaplace(0,1,0);
            nodes_[no2] -> setConstrainsLaplace(0,1,0);
            nodes_[no3] -> setConstrainsLaplace(0,1,0);
        };
        if ((boundary_[ibound] -> getConstrain(1) != 2)){
            nodes_[no1] -> setConstrainsLaplace(1,1,0);
            nodes_[no2] -> setConstrainsLaplace(1,1,0);
            nodes_[no3] -> setConstrainsLaplace(1,1,0);
        };
        
        if ((boundary_[ibound] -> getConstrain(0) == 1) || (boundary_[ibound] -> getConstrain(0) == 3)){

            //Desfazer primeira parte do if para voltar a cond. cont. constante
            // if (boundary_[ibound] -> getConstrainValue(0) <= 1.){
                
            //     typename Node::VecLocD x;
            //     x = nodes_[no1]->getCoordinates();                 
                
            //     nodes_[no1] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no2]->getCoordinates();                 
                
            //     nodes_[no2] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
                
            //     x = nodes_[no3]->getCoordinates();                 
                
            //     nodes_[no3] -> setConstrains(0,boundary_[ibound] -> 
            //                                  getConstrain(0),
            //                                  x(1) * boundary_[ibound] -> 
            //                                  getConstrainValue(0));
            //     //ate aqui
            // } else {
            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                     boundary_[ibound] -> getConstrainValue(0));
             // };
        };

        if((boundary_[ibound] -> getConstrain(1) == 1) || (boundary_[ibound] -> getConstrain(1) == 3)){
            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                     boundary_[ibound] -> getConstrainValue(1));
        };     
    };



    // For Cavity flow only
    // for (int k = 0; k < numNodes; ++k)
    // {
    //     typename Node::VecLocD x;
    //     x = nodes_[k] -> getCoordinates();
    //     if ((x(0) < 0.001) || (x(0) > 0.999))
    //     {
    //         nodes_[k] -> setConstrains(0,1.,0.);
    //         nodes_[k] -> setConstrains(1,1.,0.);
    //     }
    // }
























        //Print nodal constrains
    for (int i=0; i<numNodes; i++){

        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0) \
                  << " " << nodes_[i] -> getConstrainValue(0) \
                  << " " << nodes_[i] -> getConstrains(1) \
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 




    //Sets fluid elements and sides on interface boundaries
    for (int i=0; i<numBoundElems; i++){
        if ((boundary_[i] -> getConstrain(0) > 0) ||
            (boundary_[i] -> getConstrain(1) > 0)) {

            Boundaries::BoundConnect connectB;
            connectB = boundary_[i] -> getBoundaryConnectivity();

            for (int j=0; j<numElem; j++){
                typename Elements::Connectivity connect;
                connect = elements_[j] -> getConnectivity();
                
                int flag = 0;
                int side[3];
                for (int k=0; k<6; k++){
                    if ((connectB(0) == connect(k)) || 
                        (connectB(1) == connect(k)) ||
                        (connectB(2) == connect(k))){
                        side[flag] = k;
                        flag++;
                    };
                };
                if (flag == 3){
                    boundary_[i] -> setElement(elements_[j] -> getIndex());
                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundary_[i] -> setElementSide(0);
                        elements_[j] -> setElemSideInBoundary(0);
                    };
                    

                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundary_[i] -> setElementSide(1);
                        elements_[j] -> setElemSideInBoundary(1);
                    };
                    
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundary_[i] -> setElementSide(2);
                        elements_[j] -> setElemSideInBoundary(2);
                    };
                };
            };
        };
    };


    //Sets Viscosity, density, time step, time integration 
    //scheme and field forces values
    for (int i=0; i<numElem; i++){
        elements_[i] -> setViscosity(viscInf);
        elements_[i] -> setDensity(rhoInf);
        elements_[i] -> setTimeStep(dTime);
        elements_[i] -> setTimeIntegrationScheme(integScheme);
        elements_[i] -> setFieldForce(fieldForces);
    };

    //Closing the file
    file.close();
    if (deleteFiles)
        system((remove2 + inputFile).c_str());

    // Fluid-structure interaction data reading
    std::ifstream input(structuralInputFile);

    // Center coordinates
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> coords = splitdouble(line, " ");
    // Center initial velocity
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> vel = splitdouble(line, " ");
    // Center initial acceleration
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> accel = splitdouble(line, " ");

    Node3DOF center(0,coords,vel,accel);
    
    // Identifying nodes from the moving boundary data
    int nnodes = 0;
    std::vector<StructuralNode> FSINodes;
    for(auto nodeNumber: FSIBoundaryNodesNumbers){
        std::vector<double> coords{nodes_[nodeNumber] -> getCoordinateValue(0), nodes_[nodeNumber] -> getCoordinateValue(1)};
        FSINodes.push_back(StructuralNode(nodeNumber, coords));
        nnodes += 1;
    }

    // Constraints (1 for true, 0 for false)
    std::getline(input, line);
    std::getline(input, line);
    std::vector<bool> constraints = splitbool(line, " ");
    // Springs
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> springs = splitdouble(line, " ");
    // Masses
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> masses = splitdouble(line, " ");
    // Dampers
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> dampers = splitdouble(line, " ");
    // beta, gamma
    std::getline(input, line);
    std::getline(input, line);
    std::vector<double> aux = splitdouble(line, " ");
    // Time step at which the structure is allowed to move
    std::getline(input, line);
    std::getline(input, line);
    initialStructuralTimeStep_ = stoi(line);
    // Time step before the structure moves
    std::getline(input, line);
    std::getline(input, line);
    firstTimeStep_ = stod(line);

    structure_ = StructuralDomain(center,FSINodes,constraints,springs,masses,dampers,dTime,aux[0],aux[1]);

    return;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::readInitialValues(const std::string& inputVel,const std::string& inputPres) {

    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::ifstream inputPressure(inputPres.c_str());

    std::ifstream inputVelocity(inputVel.c_str());

    std::string line;

    for (int i = 0; i < numNodes; ++i)
    {
        double u_[2];
        double uz;
        inputVelocity >> u_[0] >> u_[1] >> uz;
        //getline(inputVelocity,line);
        nodes_[i] -> setVelocity(u_);
        //if (rank==0)std::cout << "asdasd " << i << " " << u_[0] << " " << u_[1] << " " << uz << std::endl;
    }

    for (int i = 0; i < numNodes; ++i)
    {
        double a_[2];
        double p_;
        inputPressure >> a_[0] >> a_[1] >> p_;
        //getline(inputPressure,line);
        nodes_[i] -> setPressure(p_);
        //if (rank==0)std::cout << "pressure " << i << " " << a_[0] << " " << a_[1] << " " << p_ << std::endl;
    }
    

    //  std::ifstream inputData(inputFile.c_str());
    // std::ofstream mirrorData(mirror.c_str());
    // std::ifstream file(inputMesh);
    // std::string line;
    // std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);
  


    
    

    return;
}


//------------------------------------------------------------------------------
//-------------------------SOLVE STEADY LAPLACE PROBLEM-------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveSteadyLaplaceProblem(int iterNumber, double tolerance) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                            2*numNodes, 2*numNodes,
                            50,NULL,50,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

        for (int jel = 0; jel < numElem; jel++){   
            
            //if (part_elem[jel] == rank) {
            
            //Compute Element matrix
            elements_[jel] -> getSteadyLaplace();
                            
            typename Elements::LocalMatrix Ajac;
            typename Elements::LocalVector Rhs;
            typename Elements::Connectivity connec;
            
            //Gets element connectivity, jacobian and rhs 
            connec = elements_[jel] -> getConnectivity();
            Ajac = elements_[jel] -> getJacNRMatrix();
            Rhs = elements_[jel] -> getRhsVector();
   
            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<6; i++){
                for (int j=0; j<6; j++){

                    int dof_i = 2*connec(i);
                    int dof_j = 2*connec(j);
                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,            \
                                        &Ajac(2*i  ,2*j  ),ADD_VALUES);
                    
                    dof_i = 2*connec(i)+1;
                    dof_j = 2*connec(j);
                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,            \
                                        &Ajac(2*i+1,2*j  ),ADD_VALUES);
                    
                    dof_i = 2*connec(i);
                    dof_j = 2*connec(j)+1;
                    ierr = MatSetValues(A,1,&dof_i,1,&dof_j,            \
                                        &Ajac(2*i  ,2*j+1),ADD_VALUES);
                    
                    dof_i = 2*connec(i)+1;
                    dof_j = 2*connec(j)+1;
                        ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                            &Ajac(2*i+1,2*j+1),ADD_VALUES);
                };
                                    
                //Rhs vector
            // if (fabs(Rhs(2*i  )) >= 1.e-8){
                int dofv_i = 2*connec(i);
                ierr = VecSetValues(b,1,&dofv_i,&Rhs(2*i  ),ADD_VALUES);
                //  };
                
                //if (fabs(Rhs(2*i+1)) >= 1.e-8){
                dofv_i = 2*connec(i)+1;
                ierr = VecSetValues(b,1,&dofv_i,&Rhs(2*i+1),ADD_VALUES);
                //};
            };
        };
   
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
#if defined(PETSC_HAVE_MUMPS)
        ierr = KSPSetType(ksp,KSPPREONLY);
        ierr = KSPGetPC(ksp,&pc);
        ierr = PCSetType(pc, PCLU);
#endif
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);
        
        
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        //std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            nodes_[i] -> incrementCoordinate(0,u_[0]);
            nodes_[i] -> incrementCoordinate(1,u_[1]);
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "MESH MOVING - ERROR = " << val 
                      << std::scientific <<  std::endl;
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);

        if(val <= tolerance){
            break;            
        };             
    };
    
    // for (int i=0; i<numElem; i++){
    //     elements_[i] -> computeNodalGradient();            
    // };

    if (rank == 0) {
        //Computing velocity divergent
        //      printResults(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveTransientProblem(int iterNumber, double tolerance) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //IS             rowperm       = NULL,colperm = NULL;
    //    MatNullSpace      nullsp;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::ofstream dragLift;
    std::ofstream centerDisplacements;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    centerDisplacements.open("Center_Displacements_Output.txt", std::ofstream::out | std::ofstream::app);
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };    
        
    iTimeStep = 0;



    double dTimeAux=dTime;
    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){

        for (int i = 0; i < numElem; i++){
            elements_[i] -> getParameterSUPG();
        };
    
        //Start the analysis with first order time integration and then change to the user defined
        double one = 1.;
        if (iTimeStep == 0 && (integScheme < 1. || dTimeAux != firstTimeStep_)){
            for (int i = 0; i < numElem; i++){
                elements_[i] -> setTimeIntegrationScheme(one);
                dTime = firstTimeStep_; //passo de tempo inicial
                elements_[i] -> setTimeStep(dTime);
            };
        };
        if (iTimeStep == initialStructuralTimeStep_) {
            dTime = dTimeAux; //passo de tempo permanente
            for (int i = 0; i < numElem; i++){
                elements_[i] -> setTimeStep(dTime);
            };
        };
        if (iTimeStep == initialStructuralTimeStep_ + 10) {
            for (int i = 0; i < numElem; i++){
                elements_[i] -> setTimeIntegrationScheme(integScheme);
            };
        };

        //set different  iterationNumbers 
        int iterNumber2=iterNumber;
        if(iTimeStep < 4)iterNumber2 = 10;


        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        for (int i = 0; i < numNodes; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodes_[i] -> getVelocity(0);
            u[1] = nodes_[i] -> getVelocity(1);
            
            uprev[0] = nodes_[i] -> getPreviousVelocity(0);
            uprev[1] = nodes_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodes_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodes_[i] -> setPreviousVelocity(u);

        };

        double duNorm=100.;
         
        for (int inewton = 0; inewton < iterNumber2; inewton++){
            boost::posix_time::ptime t1 =                             
                               boost::posix_time::microsec_clock::local_time();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2*numNodes+numNodes, 2*numNodes+numNodes,
                                100,NULL,300,NULL,&A); 
            CHKERRQ(ierr);
            
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
            
            //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;

            for (int jel = 0; jel < numElem; jel++){   
                
                //if (part_elem[jel] == rank) {
                    //Compute Element matrix
                    elements_[jel] -> getTransientNavierStokes();

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elements_[jel] -> getConnectivity();
                    Ajac = elements_[jel] -> getJacNRMatrix();
                    Rhs = elements_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    //Matrix K and C
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                        
                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodes + connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };
                        
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i  ),
                                                ADD_VALUES);
                        };
                        
                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i)+1;
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodes + connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                //};
            }; //Elements
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            

            // Mat Aperm;
            // MatGetOrdering(A,MATORDERINGRCM,&rowperm,&colperm);
            // MatPermute(A,rowperm,colperm,&Aperm);
            // VecPermute(b,colperm,PETSC_FALSE);
            // MatDestroy(&A);
            // A    = Aperm;    

            //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //MatView(A,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
            //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            
            // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
            //                         500);CHKERRQ(ierr);
            
            // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            
            // // ierr = KSPGetPC(ksp,&pc);
            
            // // ierr = PCSetType(pc,PCNONE);
            
            // // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);

            // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);
            
            // //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            

        // //   //   ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
        // // // ierr = MatSetNullSpace(A, nullsp);
        // // // ierr = MatNullSpaceDestroy(&nullsp);

   

#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
#endif          
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);



            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

            ierr = KSPGetTotalIterations(ksp, &iterations);            

            // VecPermute(u,rowperm,PETSC_TRUE);

            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double p_;
            duNorm = 0.;
            double dpNorm = 0.;
            Ione = 1;
            
            for (int i = 0; i < numNodes; ++i){
                //if (nodes_[i] -> getConstrains(0) == 0){
                    Ii = 2*i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(0,val);
                    //}; 
                
                    //if (nodes_[i] -> getConstrains(1) == 0){
                    Ii = 2*i+1;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(1,val);
                    //};
            };
            
            for (int i = 0; i<numNodes; i++){
                Ii = 2*numNodes+i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                dpNorm += val*val;
                nodes_[i] -> incrementPressure(p_);
            };
            
            //Computes the solution vector norm
            //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
   
            boost::posix_time::ptime t2 =                               \
                               boost::posix_time::microsec_clock::local_time();
         
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;

                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                          << " " << sqrt(dpNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds()/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);

            if ((iTimeStep>2)&&(sqrt(duNorm) <= tolerance)) {
                break;
            };

            
        };//Newton-Raphson

        // Compute and print drag and lift coefficients
        if (computeDragAndLift){
            dragAndLiftCoefficients(dragLift);
        };
        

        if (iTimeStep >= initialStructuralTimeStep_){
            std::vector<double> F(3,0.0);
            computeFSIForces(F);
            if (rank == 0){
                std::cout << "FORCES ACTING ON THE STRUCTURE: " << F[0] << " " << F[1] << " " << F[2] << std::endl;
            }

            structure_.updateCenterPosition(F);
            if (rank == 0){
                std::cout << "STRUCTURE'S CENTER POSITION: " << structure_.getCenter().getCurrentPosition(0) << " " << structure_.getCenter().getCurrentPosition(1) << " " << structure_.getCenter().getCurrentPosition(2) << std::endl; 
            }
            structure_.updateBoundary();
        }

        //Printing results
        printResults(iTimeStep);
        printFSIBoundary(iTimeStep,"FSI_Boundary_Output");
        printStructureCenterDisplacements(iTimeStep,centerDisplacements);
        
    };

    centerDisplacements.close();
    dragLift.close();
    
    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveTransientProblemMoving(int iterNumber, double tolerance) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    //    MatNullSpace      nullsp;
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };

    // Set element mesh moving parameters
    double vMax = 0., vMin = 1.e10;
    for (int i = 0; i < numElem; i++){
        double v = elements_[i] -> getJacobian();
        if (v > vMax) vMax = v;
        if (v < vMin) vMin = v;
    };
    for (int i = 0; i < numElem; i++){
        double v = elements_[i] -> getJacobian();
        double eta = 1 + (1. - vMin / vMax) / (v / vMax);
        elements_[i] -> setMeshMovingParameter(eta);
    };

    iTimeStep = 0;
    double dTimeAux = dTime;

    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){

        for (int i = 0; i < numElem; i++){
            elements_[i] -> getParameterSUPG();
        };

        //Start the analysis with first order time integration and then change to the user defined
        if (integScheme < 1.){
            double one = 1.;
            if (iTimeStep == 0){
                for (int i = 0; i < numElem; i++){
                    elements_[i] -> setTimeIntegrationScheme(one);
                    dTime = firstTimeStep_; //passo de tempo inicial
                    elements_[i] -> setTimeStep(dTime);
                };
            };
            if (iTimeStep == initialStructuralTimeStep_) {
                for (int i = 0; i < numElem; i++){
                    elements_[i] -> setTimeIntegrationScheme(integScheme);
                    dTime = dTimeAux; //passo de tempo permanente
                    elements_[i] -> setTimeStep(dTime);
                };
            };
        };

        //set different  iterationNumbers 
        int iterNumber2=iterNumber;
        if(iTimeStep < 4)iterNumber2 = 10;
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        for (int i = 0; i < numNodes; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodes_[i] -> getVelocity(0);
            u[1] = nodes_[i] -> getVelocity(1);
            
            uprev[0] = nodes_[i] -> getPreviousVelocity(0);
            uprev[1] = nodes_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodes_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodes_[i] -> setPreviousVelocity(u);
        };


        // Moving boundary
        for (int i = 0; i < numNodes; i++){
            typename Node::VecLocD x;
            
            x = nodes_[i] -> getCoordinates();
            nodes_[i] -> setPreviousCoordinates(0,x(0));
            nodes_[i] -> setPreviousCoordinates(1,x(1));
        };

        for (int i=0; i < numBoundElems; i++){
            if (boundary_[i] -> getConstrain(0) == 3){
                //std::cout << "asasa " << i << std::endl;
                Boundaries::BoundConnect connectB;
                connectB = boundary_[i] -> getBoundaryConnectivity();
                int no1 = connectB(0);
                int no2 = connectB(1);
                int no3 = connectB(2);
                
                typename Node::VecLocD x;
                x = nodes_[no1]->getCoordinates();
                x(1) -= .1 * dTime*0;
                nodes_[no1] -> setUpdatedCoordinates(x);

                x = nodes_[no2]->getCoordinates();
                x(1) -= .1 * dTime*0;
                nodes_[no2] -> setUpdatedCoordinates(x);

                x = nodes_[no3]->getCoordinates();
                x(1) -= .1 * dTime*0;
                nodes_[no3] -> setUpdatedCoordinates(x);
            };
        };

        solveSteadyLaplaceProblem(1, 1.e-6);

        for (int i=0; i< numNodes; i++){
            typename Node::VecLocD x,xp;
            double u[2];
            
            x = nodes_[i] -> getCoordinates();
            xp = nodes_[i] -> getPreviousCoordinates();
            
            u[0] = (x(0) - xp(0)) / dTime;
            u[1] = (x(1) - xp(1)) / dTime;

            nodes_[i] -> setMeshVelocity(u);
        };


        double duNorm=100.;
        
        for (int inewton = 0; inewton < iterNumber2; inewton++){
            boost::posix_time::ptime t1 =                             
                               boost::posix_time::microsec_clock::local_time();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                2*numNodes+numNodes, 2*numNodes+numNodes,
                                100,NULL,300,NULL,&A); 
            CHKERRQ(ierr);
            
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
            
            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodes);
            CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
            
            for (int jel = 0; jel < numElem; jel++){   
                
                //if (part_elem[jel] == rank) {
                    //Compute Element matrix
                    elements_[jel] -> getTransientNavierStokes();

                    typename Elements::LocalMatrix Ajac;
                    typename Elements::LocalVector Rhs;
                    typename Elements::Connectivity connec;
                    
                    //Gets element connectivity, jacobian and rhs 
                    connec = elements_[jel] -> getConnectivity();
                    Ajac = elements_[jel] -> getJacNRMatrix();
                    Rhs = elements_[jel] -> getRhsVector();
                    
                    //Disperse local contributions into the global matrix
                    //Matrix K and C
                    for (int i=0; i<6; i++){
                        for (int j=0; j<6; j++){
                            if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                    &Ajac(2*i  ,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,2*j+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * connec(j) + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,2*j+1),
                                                    ADD_VALUES);
                            };
                        
                            //Matrix Q and Qt
                            if (fabs(Ajac(2*i  ,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i  ,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i  )) >= 1.e-15){
                                int dof_i = 2 * connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i  ),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(2*i+1,12+j)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(2*i+1,12+j),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+j,2*i+1)) >= 1.e-15){
                                int dof_i = 2 * connec(i) + 1;
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                                    &Ajac(12+j,2*i+1),
                                                    ADD_VALUES);
                            };
                            if (fabs(Ajac(12+i,12+j)) >= 1.e-15){
                                int dof_i = 2 * numNodes + connec(i);
                                int dof_j = 2 * numNodes + connec(j);
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &Ajac(12+i,12+j),
                                                    ADD_VALUES);
                            };
                        };
                        
                        //Rhs vector
                        if (fabs(Rhs(2*i  )) >= 1.e-15){
                            int dof_i = 2 * connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i  ),
                                                ADD_VALUES);
                        };
                        
                        if (fabs(Rhs(2*i+1)) >= 1.e-15){
                            int dof_i = 2 * connec(i)+1;
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(2*i+1),
                                                ADD_VALUES);
                        };
                        if (fabs(Rhs(12+i)) >= 1.e-15){
                            int dof_i = 2 * numNodes + connec(i);
                            ierr = VecSetValues(b, 1, &dof_i, &Rhs(12+i),
                                                ADD_VALUES);
                        };
                    };
                //};
            }; //Elements
            
            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
            
            // MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            // ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            



            // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
            //                         500);CHKERRQ(ierr);
            
            //ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            
            // ierr = KSPGetPC(ksp,&pc);
            
            // ierr = PCSetType(pc,PCJACOBI);
            
            // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);

            //ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);
            
            // //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
            

            // ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&nullsp);
            // ierr = MatSetNullSpace(A, nullsp);
            // ierr = MatNullSpaceDestroy(&nullsp);

 

#if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
            //      MatMumpsSetIcntl(A,25,-1);
#endif          
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);


            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);

            ierr = KSPGetTotalIterations(ksp, &iterations);            

            //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
            
            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);
            CHKERRQ(ierr);
            
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            double p_;
            duNorm = 0.;
            double dpNorm = 0.;
            Ione = 1;
            
            for (int i = 0; i < numNodes; ++i){
                //if (nodes_[i] -> getConstrains(0) == 0){
                    Ii = 2*i;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(0,val);
                    //}; 
                
                    //if (nodes_[i] -> getConstrains(1) == 0){
                    Ii = 2*i+1;
                    ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                    duNorm += val*val;
                    nodes_[i] -> incrementVelocity(1,val);
                    //};
            };
            
            for (int i = 0; i<numNodes; i++){
                Ii = 2*numNodes+i;
                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
                p_ = val;
                dpNorm += val*val;
                nodes_[i] -> incrementPressure(p_);
            };
            
            //Computes the solution vector norm
            //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
   
            boost::posix_time::ptime t2 =                               \
                               boost::posix_time::microsec_clock::local_time();
         
            if(rank == 0){
                boost::posix_time::time_duration diff = t2 - t1;

                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(duNorm) 
                          << " " << sqrt(dpNorm)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds()/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);

            if ((iTimeStep>2)&&(sqrt(duNorm) <= tolerance)) {
                break;
            };

            //Updates SUPG Parameter
            // for (int i = 0; i < numElem; i++){
            //     elements_[i] -> getParameterSUPG();
            // };
            
        };//Newton-Raphson

        // Compute and print drag and lift coefficients
        if (computeDragAndLift){
            dragAndLiftCoefficients(dragLift);
        };

        if (printVorticity){
            for (int i = 0; i < numNodes; i++){
                nodes_[i] -> clearVorticity();
            };
            for (int jel = 0; jel < numElem; jel++){
                elements_[jel] -> computeVorticity();
            };

            for (int i = 0; i < numNodes; ++i){
                double vort = nodes_[i] -> getVorticity();
                double signal = vort / fabs(vort);
                int root;
                struct { 
                    double val; 
                    int   rank; 
                } in, out; 

                in.val = fabs(vort);
                in.rank = rank;

                MPI_Reduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,root,PETSC_COMM_WORLD);
                MPI_Bcast(&out.val,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
                MPI_Bcast(&out.rank,1,MPI_INT,0,PETSC_COMM_WORLD);
                MPI_Bcast(&signal,1,MPI_DOUBLE,out.rank,PETSC_COMM_WORLD);

                vort = out.val * signal;

                nodes_[i] -> setVorticity(vort); 
            }


        };

        //Printing results
        printResults(iTimeStep);
        
    };

   
    return 0;
};



#endif
