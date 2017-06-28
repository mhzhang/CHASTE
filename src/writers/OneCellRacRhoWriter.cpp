
/*

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#include "OneCellRacRhoWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VertexElement.hpp"
#include <VertexBasedCellPopulation.hpp>
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::OneCellRacRhoWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("OneCellData.csv")
{
}

// Writing Header
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	
	*this->mpOutStream << "# TimeStamp,Cell_ID,Cell_Area,Target_Area,Cell_Perimeter,num_neighbours,num_edges,rac_concentration,rho_concentration\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << ",";
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellRacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		double cell_id = cell_iter->GetCellId();
		if (cell_id == 1){
			double volume = cell_iter->GetCellData()->GetItem("volume");
			unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
			double perimeter = pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(elem_index);
			std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringLocationIndices(*cell_iter);
			int num_neighbours = neighbour_indices.size();
			VertexElement < SPACE_DIM, SPACE_DIM > *VertexElement = pCellPopulation->GetElementCorrespondingToCell(*cell_iter);
			int num_edges = VertexElement->GetNumNodes();
			double rac = cell_iter->GetCellData()->GetItem("F");
			double rho = cell_iter->GetCellData()->GetItem("W");
			double target_area = cell_iter->GetCellData()->GetItem("target area");
			//double area = cell_iter->GetCellData()->GetItem("AREA");
		    *this->mpOutStream  << ","<< cell_id;
			*this->mpOutStream  << ","<< volume;
			*this->mpOutStream  << ","<< target_area;
			*this->mpOutStream  << ","<< perimeter;
			*this->mpOutStream  << ","<< num_neighbours;
			*this->mpOutStream  << ","<< num_edges;
			*this->mpOutStream  << ","<< rac;
			*this->mpOutStream  << ","<< rho;
			//*this->mpOutStream  << ","<< area;
			
		}

		}
}

// Explicit instantiation
template class OneCellRacRhoWriter<1,1>;
template class OneCellRacRhoWriter<1,2>;
template class OneCellRacRhoWriter<2,2>;
template class OneCellRacRhoWriter<1,3>;
template class OneCellRacRhoWriter<2,3>;
template class OneCellRacRhoWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OneCellRacRhoWriter)
