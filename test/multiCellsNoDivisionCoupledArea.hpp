/*
TEMPLATE FOR GTPASE SIMULATIONS WITH COUPLED CELL AREA

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
//#include "TargetAreaModifierCoupledArea.hpp"
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
#include "ODEParameterAreaModifier.hpp"
#include "ShapeWriter.hpp"
#include "CsvWriter.hpp"
#include "NumNeighboursWriter.hpp"
#include "XMLCellWriter.hpp"
//#include "RacRhoWriter.hpp"
#include "OneCellGTPaseWriter.hpp"

#include "ODESRNCoupledArea.hpp"
#include "AbstractOdeSrnModel.hpp"

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
//#include "AbstractOdeSrnModel.hpp"

#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include <cmath>


class multiCellsNoDivisionCoupled : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedMonolayer() throw (Exception)
    {

        HoneycombVertexMeshGenerator generator(20, 20); //20x20 = 400cells
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
			double G_conc = (randGenerator->ranf());

			std::vector<double> initial_conditions;

			initial_conditions.push_back(G_conc); //G, used 1 previously, now set to random
			initial_conditions.push_back(0.8); //target_area
			initial_conditions.push_back(0.866025); //area
			p_srn_model->SetInitialConditions(initial_conditions);
			
			
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-(double)i - 2.0); // So all out of M phase
            p_cycle_model->SetQuiescentVolumeFraction(1.0);
            p_cycle_model->SetEquilibriumVolume(1.0);

			

            CellPtr p_cell(new Cell(p_state, p_cycle_model, p_srn_model));

            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();


			
			
            cells.push_back(p_cell);
			counter++;
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
		//cell_population.AddCellWriter<XMLCellWriter>();
		cell_population.AddPopulationWriter<CsvWriter>();
		cell_population.AddPopulationWriter<ShapeWriter>();
		cell_population.AddPopulationWriter<NumNeighboursWriter>();
		//cell_population.AddPopulationWriter<RacRhoWriter>();
		cell_population.AddPopulationWriter<OneCellGTPaseWriter>();
		cell_population.AddCellWriter<XMLCellWriter>();
		
		
		//Following code used to add cell label, not working as of now.
		/*
		boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cell_iter->AddCellProperty(p_label);
            }
        }
		*/
		

        OffLatticeSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("testDifferentialAdhesionHighFinal"); //output directory name
        simulator.SetSamplingTimestepMultiple(200);//200
		simulator.SetDt(0.01); //0.01 for 50x50
        simulator.SetEndTime(120.0); //2500 for proper simulation length
		
        
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
		
		MAKE_PTR(ODEParameterAreaModifier<2>, p_ODE_modifier);
		simulator.AddSimulationModifier(p_ODE_modifier);
		

		
		//For low adhesion, want to minimize difference between CellBoundaryAdhesion and CellCellAdhesion
		//CellBoundaryAdhesion set to 1 by default
		//CellCellAdhesion set to 0.5 by default
        //MAKE_PTR(NagaiHondaForce<2>, p_force);
		MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
		//Following code used for differential adhesion, not working as of now.
		
		p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
		p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12); //low Boundary Adhesion
		p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(40.0);
		p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1); //high adhesion between 2 non labelled cells
		p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(6); //low adhesion between labelled and non labelled cell
		p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(3); //normal adhesion between 2 labelled cells
        simulator.AddForce(p_force);

        

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();


    }


};

#endif /*MULTICELLSNODIVISIONCOUPLED_HPP_*/
