/*
COUPLED RAC AND RHO SIMULATION

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/
#ifndef MULTICELLSNODIVISIONCOUPLED_HPP_
#define MULTICELLSNODIVISIONCOUPLED_HPP_


#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "VoronoiDataWriter.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"

//#include "FakePetscSetup.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimulationTime.hpp"
#include "CellLabel.hpp"
#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "ContactInhibitionCellCycleModel.hpp"

#include "VolumeTrackingModifier.hpp"
#include "ODEParameterAverageModifier.hpp"
#include "ShapeWriter.hpp"
#include "CsvWriter.hpp"
#include "NumNeighboursWriter.hpp"
#include "RacRhoWriter.hpp"
#include "OneCellRacRhoWriter.hpp"

#include "ODESRNCoupledAverage.hpp"
#include "AbstractOdeSrnModel.hpp"

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractOdeSrnModel.hpp"

#include <cmath>

class multiCellsNoDivisionCoupled : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Test 1 - a basic vertex-based simulation ==
     *
     * In the first test, we run a simple vertex-based simulation of an epithelial monolayer.
     * Each cell in the simulation is assigned a simple stochastic cell-cycle model, the cells will divide randomly and never stop proliferating.
     */
    void TestVertexBasedMonolayer() throw (Exception)
    {
        /* The first thing we define is a 2D (specified by the <2,2>) mesh which holds the spatial information of the simulation. To do this we use one of a
         * number of {{{MeshGenerators}}}.*/
        HoneycombVertexMeshGenerator generator(20, 20);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
		RandomNumberGenerator* randGenerator = RandomNumberGenerator::Instance();
		randGenerator->Reseed(1);
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
			ODESrnModel* p_srn_model = new ODESrnModel;
			
			std::vector<double> initial_conditions;

			double rac_conc = (randGenerator->ranf());
			double rho_conc = (randGenerator->ranf());
			initial_conditions.push_back(rac_conc); //rac
			initial_conditions.push_back(rho_conc); //rho
			initial_conditions.push_back(rac_conc); //rac_average
			initial_conditions.push_back(rho_conc); //rho_average
			p_srn_model->SetInitialConditions(initial_conditions);
			
			
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-(double)i - 2.0); // So all out of M phase
            p_cycle_model->SetQuiescentVolumeFraction(1.0);
            p_cycle_model->SetEquilibriumVolume(1.0);

			

            CellPtr p_cell(new Cell(p_state, p_cycle_model, p_srn_model));
			//p_srn_model->ResetForDivision();
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
			p_cell->GetCellData()->SetItem("F", rac_conc);
			p_cell->GetCellData()->SetItem("W", rho_conc);
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
		//cell_population.AddCellWriter<XMLCellWriter>();
		cell_population.AddPopulationWriter<CsvWriter>();
		cell_population.AddPopulationWriter<ShapeWriter>();
		cell_population.AddPopulationWriter<NumNeighboursWriter>();
		cell_population.AddPopulationWriter<RacRhoWriter>();
		cell_population.AddPopulationWriter<OneCellRacRhoWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
		//simulator.AddCellPopulationBoundaryCondition(p_bc);
        simulator.SetOutputDirectory("3x3run_no_division_rac1rho1_coupled_average_random_init_2");
        simulator.SetSamplingTimestepMultiple(10);
		simulator.SetDt(0.01);
        simulator.SetEndTime(80.0); //120
		
        
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
		
		MAKE_PTR(ODEParameterAverageModifier<2>, p_ODE_modifier);
		simulator.AddSimulationModifier(p_ODE_modifier);
		
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        /* In order to specify how cells move around we create a "shared pointer" to a
         * {{{Force}}} object and pass it to the {{{OffLatticeSimulation}}}. This is done using the MAKE_PTR macro as follows.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* A {{{NagaiHondaForce}}} has to be used together with a child class of {{{AbstractTargetAreaModifier}}}.
         * This modifies the target area of individual cells and thus alters the relative forces
         * between neighbouring cells.
         */
        

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();



    }


};

#endif /*MULTICELLSNODIVISIONCOUPLED_HPP_*/
