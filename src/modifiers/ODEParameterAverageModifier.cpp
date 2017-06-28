/*

A MODIFIER TO COUPLE THE NEIGHBOURING CELLS' AVERAGE RAC AND RHO CONCNTRATIONS
Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#include "ODEParameterAverageModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractOdeSrnModel.hpp"
#include "AbstractSrnModel.hpp"

template<unsigned DIM>
ODEParameterAverageModifier<DIM>::ODEParameterAverageModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ODEParameterAverageModifier<DIM>::~ODEParameterAverageModifier()
{
}

template<unsigned DIM>
void ODEParameterAverageModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAverageModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAverageModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
	//VertexBasedCellPopulation<DIM>* pCellPopulation = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
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
		double rac_conc = 0;
		double rho_conc = 0;
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
		for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
		neighbour_iter != neighbour_indices.end();
		++neighbour_iter){
			unsigned neighbour_index = *neighbour_iter;
			CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);
			rac_conc += p_neighbour_cell->GetCellData()->GetItem("F");
			rho_conc += p_neighbour_cell->GetCellData()->GetItem("W");
		}
		int num_neighbours = neighbour_indices.size();
		double avg_rac_conc = rac_conc/num_neighbours;
		double avg_rho_conc = rho_conc/num_neighbours;
        AbstractSrnModel* model_ptr_base = cell_iter->GetSrnModel();
		AbstractOdeSrnModel* model_ptr = dynamic_cast<AbstractOdeSrnModel*>(model_ptr_base);
		model_ptr->GetStateVariables()[2]=avg_rac_conc;
		model_ptr->GetStateVariables()[3]=avg_rho_conc;
    }
}

template<unsigned DIM>
void ODEParameterAverageModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ODEParameterAverageModifier<1>;
template class ODEParameterAverageModifier<2>;
template class ODEParameterAverageModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ODEParameterAverageModifier)

