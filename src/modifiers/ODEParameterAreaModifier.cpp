/*

A MODIFIER TO COUPLE THE CELL'S AREA TO THE ODESRN SYSTEM
Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#include "ODEParameterAreaModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractOdeSrnModel.hpp"
#include "AbstractSrnModel.hpp"

template<unsigned DIM>
ODEParameterAreaModifier<DIM>::ODEParameterAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ODEParameterAreaModifier<DIM>::~ODEParameterAreaModifier()
{
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    /**
     * This hack is needed because in the case of a MeshBasedCellPopulation in which
     * multiple cell divisions have occurred over one time step, the Voronoi tessellation
     * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
     * tessellation here, an assertion may trip as we try to access a Voronoi element
     * whose index exceeds the number of elements in the out-of-date tessellation.
     *
     * \todo work out how to properly fix this (#1986)
     */
    if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the volume of this cell
        //double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);
		double cell_volume = cell_iter->GetCellData()->GetItem("volume");

		
        AbstractSrnModel* model_ptr_base = cell_iter->GetSrnModel();
		AbstractOdeSrnModel* model_ptr = dynamic_cast<AbstractOdeSrnModel*>(model_ptr_base);
		model_ptr->GetStateVariables()[2]=cell_volume;
    }
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ODEParameterAreaModifier<1>;
template class ODEParameterAreaModifier<2>;
template class ODEParameterAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ODEParameterAreaModifier)

