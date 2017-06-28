/*
AN ODESRN SYSTEM WITHOUT ANY COUPLING
Author: MoHan Zhang
Contact: mohan_z@hotmail.com

 */

#ifndef ODESRN_HPP_
#define ODESRN_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header includes the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/* The next header defines a base class for ode-based SRN models.
 * Our new SRN model will inherit from this abstract class. */
#include "AbstractOdeSrnModel.hpp"

/* These headers specify the methods to solve the ODE system.*/
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/* This header specifies the ODE solvers.*/
#include "CellCycleModelOdeSolver.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based Chaste
 * tutorials. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#include <cmath>
/*
 * == Defining the SRN model and ODE system classes ==
 *
 * As an example, let us consider a SRN model in which we solve a simple ODE
 * dx/dt = -0.25*y
 * dy/dt = x
 * This has exact solution x = A cos 0.5t + B sin 0.5t
 * where A and B are determined by the initial condions.
 *
 * To implement this model we define a new SRN model, {{{MySrnModel}}},
 * which inherits from {{{AbstractOdeSrnModel}}} and
 * contains a {{{MyOdeSystem}}}.
 *
 * Note that usually this code would be separated out into a separate declaration in
 * a .hpp file and definition in a .cpp file.
 */
class ODESRN : public AbstractOdeSystem
{
private:
    friend class boost::serialization::access;
    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the ODE system (and therefore the SRN model) object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the ODE system
     * using the serialization code defined in the base class
     * {{{AbstractOdeSystem}}}.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:
    ODESRN() : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<ODESRN>::Instance();
    }

	//Variables and equations in the ODE system
    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
		double da=1;
		double db=1;
		double n=3;
		double m=3;
		double IW=3;
		double IF=3;
        rDY[0] = IF/(1+pow(rY[1],m))-db*rY[0];
        rDY[1] = IW/(1+pow(rY[0],n))-da*rY[1];
    }
};

/* As in the ODE tutorials we need to define the ODE system information.
 */
template<>
void OdeSystemInformation<ODESRN>::Initialise()
{
    this->mVariableNames.push_back("F");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(10);

    this->mVariableNames.push_back("W");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(10);

    this->mInitialised = true;
}


class ODESrnModel : public AbstractOdeSrnModel
{
private:

    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the SRN model object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the SRN
     * model using the serialization code defined in the base class
     * {{{AbstractOdeSrnModel}}}.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

    /* The first public method is a constructor, which just calls the base
     * constructor.  Note you can include an optional argument to specify the ODE solver.*/
public:

    ODESrnModel()
        : AbstractOdeSrnModel(2, boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {

        mpOdeSolver = CellCycleModelOdeSolver<ODESrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
		//Time step for ODE solver, should match that in the simulation
        SetDt(0.1);

        assert(mpOdeSolver->IsSetUp());
    }

    /* The second public method overrides {{{CreateSrnModel()}}}. This is a
     * builder method to create new copies of the SRN model. We first create
     * a new SRN model, then set each member variable of the new SRN
     * model that inherits its value from the parent.
     * Finally we call the parent class method to set the current state variables
     * as initial conditions.
     */
    AbstractSrnModel* CreateSrnModel()
    {
        ODESrnModel* p_model = new ODESrnModel();

        p_model->SetOdeSystem(new ODESRN);

        return AbstractOdeSrnModel::CreateSrnModel(p_model);
    }

    /* The third public method overrides {{{Initialise()}}}. */
    void Initialise()
    {
        AbstractOdeSrnModel::Initialise(new ODESRN);
    }

    void SimulateToCurrentTime()
    {
        // run the ODE simulation as needed
        AbstractOdeSrnModel::SimulateToCurrentTime();

        /* this line outputs the ODE system variable to {{{CellData}}}. */
        mpCell->GetCellData()->SetItem("F",mpOdeSystem->rGetStateVariables()[0]);
		mpCell->GetCellData()->SetItem("W",mpOdeSystem->rGetStateVariables()[1]);
    }
	
	void ResetForDivision(){
		AbstractOdeSrnModel::ResetForDivision();
	    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
	    for (unsigned i=0; i<2; i++)
	    {
	        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
	    }
	}

};

/* We need to include the next block of code if you want to be able to archive (save or load)
 * the SRN model object in a cell-based simulation. It is also required for writing out
 * the parameters file describing the settings for a simulation - it provides the unique
 * identifier for our new SRN model. Thus every SRN model class must provide this,
 * or you'll get errors when running simulations. */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ODESRN)
CHASTE_CLASS_EXPORT(ODESrnModel)


/* Since we're defining the new SRN model and ODEs within the test file, we need to include the
 * following stanza as well, to make the code work with newer versions of the Boost libraries.
 * Normally the above export declaration would occur in the SRN model's .hpp file, and
 * the following lines would appear in the .cpp file.  See ChasteGuides/BoostSerialization for
 * more information.
 */
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ODESRN)
CHASTE_CLASS_EXPORT(ODESrnModel)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(ODESrnModel)

#endif /*ODESRN_HPP_*/
