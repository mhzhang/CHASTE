/*

SIMPLE CONTACT INHIBITION MODEL

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#ifndef TESTCELLBASEDDEMOTUTORIAL_HPP_
#define TESTCELLBASEDDEMOTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and cell-based simulations in Chaste =
 *
 * == Introduction ==
 *
 * This tutorial is designed to give you a quick introduction to running cell-based
 * simulations in Chaste. Full details are postponed until later tutorials.
 *
 * We begin with a simple monolayer simulation and see how to:
 *   * change the cell-level model;
 *   * how to impose boundaries;
 *   * how to impose periodic conditions;
 *   * how to specify how to remove cells; and
 *   * how to change cell-cycle models.
 *
 * == The test ==
 *
 * We begin by including the necessary header files. These will be described in detail in
 * subsequent cell-based tutorials.
 */
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
#include "ShapeWriter.hpp"
#include "CsvWriter.hpp"
#include "NumNeighboursWriter.hpp"
#include "RacRhoWriter.hpp"

#include "ODESRN.hpp"
#include "AbstractOdeSrnModel.hpp"

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "AbstractOdeSrnModel.hpp"

#include <cmath>
/* Next, we define the test class which inherits from {{{AbstractCellBasedTestSuite}}}.
 * We inherit from {{{AbstractCellBasedTestSuite}}} rather than {{{CxxTest::TestSuite}}} directly because
 * this class sets up and destroys some singleton objects for us. Singletons are objects that we want to exist only
 * once in each simulation and will be covered in detail in later tutorials.
 * Since we are using {{{AbstractCellBasedTestSuite}}} the singleton {{{SimulationTime}}} is initialised to zero at the beginning of the test and destroyed at the end
 * of the test; {{{RandomNumberGenerator}}} is re-seeded with zero at the beginning and destroyed at the end of the test;
 * and {{{CellPropertyRegistry}}} (which stores {{{CellProperties}}}, you learn about these in a later tutorial
 * [wiki:UserTutorials/CreatingAndUsingANewCellProperty]) is cleared at the beginning of the test.
 * This makes for cleaner code.
 */
class TestCellBasedDemoTutorial : public AbstractCellBasedTestSuite
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
        HoneycombVertexMeshGenerator generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
			ODESrnModel* p_srn_model = new ODESrnModel;
			
			std::vector<double> initial_conditions;
			//initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
			//initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
			initial_conditions.push_back(10);
			initial_conditions.push_back(10);
			p_srn_model->SetInitialConditions(initial_conditions);
			
			
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-(double)i - 2.0); // So all out of M phase
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model, p_srn_model));
			//p_srn_model->ResetForDivision();
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
		//cell_population.AddCellWriter<XMLCellWriter>();
		cell_population.AddPopulationWriter<CsvWriter>();
		cell_population.AddPopulationWriter<ShapeWriter>();
		cell_population.AddPopulationWriter<NumNeighboursWriter>();
		cell_population.AddPopulationWriter<RacRhoWriter>();
		
		/* Boundary conditions */
		//c_vector<double,2> point = zero_vector<double>(2);
        //c_vector<double,2> normal = zero_vector<double>(2);
        //normal(1) = -1.0;
		//MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
		

        /*
         * We now create an {{{OffLatticeSimulation}}} object and pass in the
         * {{{CellPopulation}}}. We also set some options on the simulation
         * like output directory, output multiple (so we don't visualize every
         * timestep), and end time.
         */
        OffLatticeSimulation<2> simulator(cell_population);
		//simulator.AddCellPopulationBoundaryCondition(p_bc);
        simulator.SetOutputDirectory("testrun5");
        simulator.SetSamplingTimestepMultiple(10);
		simulator.SetDt(0.01);
        simulator.SetEndTime(80.0); //120
		
		c_vector<double,2> point = zero_vector<double>(2);
		c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        point(0) = 6.0;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

		point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
		
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
        /* In order to specify how cells move around we create a "shared pointer" to a
         * {{{Force}}} object and pass it to the {{{OffLatticeSimulation}}}. This is done using the MAKE_PTR macro as follows.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* A {{{NagaiHondaForce}}} has to be used together with a child class of {{{AbstractTargetAreaModifier}}}.
         * This modifies the target area of individual cells and thus alters the relative forces
         * between neighbouring cells.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. If different simulation input parameters are being explored
         * the lines should be removed.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

   
    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/CellBasedDemo7/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     */
};

#endif /*TESTCELLBASEDDEMOTUTORIAL_HPP_*/
