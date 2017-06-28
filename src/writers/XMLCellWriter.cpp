
/*

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#include <boost/regex.hpp>
#include "XMLCellWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimulationTime.hpp"
//#include "CellType.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexMesh.hpp"

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::XMLCellWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cell_data.xml")
{
	this->mVtkCellDataName = "XML_dummy_attribute";
};


/*
    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::SetSimulationInfo(std::string simulationType="", std::string typeList="",
        std::string axisDivision="", std::string cellCycleModel="", std::string extraSimInfo="")
{
    mSimulationType = simulationType;
    mTypeList = typeList;
    mAxisDivision = axisDivision;
    mCellCycleModel = cellCycleModel;
    mExtraSimInfo = extraSimInfo;
}
*/


    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::OpenOutputFile(OutputFileHandler& rOutputFileHandler)
{
    AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>::OpenOutputFile(rOutputFileHandler);
	/*
    *this->mpOutStream << "<simulation simulation_type=\"" << mSimulationType << "\" ";
    *this->mpOutStream << "type_list=\"" << mTypeList << "\" ";
    *this->mpOutStream << "axis_division=\"" << mAxisDivision << "\" ";
    *this->mpOutStream << "cell_cycle_model=\"" << mCellCycleModel << "\" ";
    *this->mpOutStream << "time_step=\"" << SimulationTime::Instance()->GetTimeStep() << "\" ";
    *this->mpOutStream << "extra_sim_info=\"" << mExtraSimInfo << "\"";
    *this->mpOutStream << ">\n";
	*/
}


    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::WriteNewline()
{
    *this->mpOutStream << "</time>\n";
}

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::WriteTimeStamp()
{
    *this->mpOutStream << "<time t=\"" << SimulationTime::Instance()->GetTime() << "\" tau=\"" << SimulationTime::Instance()->GetTimeStepsElapsed() << "\">\n";
}

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	//CellID
    *this->mpOutStream << "<cell ";
	unsigned cell_id = pCell->GetCellId();
	*this->mpOutStream << "cell_id=\"" << cell_id << "\" ";
	
	//Centroids
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    *this->mpOutStream << "x=\"" << centre_location[0] << "\" ";
    *this->mpOutStream << "y=\"" << centre_location[1] << "\" ";

	//Area
	double volume = pCell->GetCellData()->GetItem("volume");
    if (volume < DBL_MAX)   // Only write cells with finite volume (avoids a case for boundary cells in MeshBasedCellPopulation)
    {
        *this->mpOutStream << "area=\"" << volume << "\" ";
    }
	
	//Label
	if (pCell->HasCellProperty<CellLabel>()){
		*this->mpOutStream << "CellLabel=\"" << 1 << "\" ";
	}
	
	//Target Area
	double target_area = pCell->GetCellData()->GetItem("target area");
	*this->mpOutStream << "target_area=\"" << target_area << "\" ";
	
	//Area from ODE
	double ODE_area = pCell->GetCellData()->GetItem("AREA");
	*this->mpOutStream << "ODE_area=\"" << ODE_area << "\" ";
	
	//GTPase Concentration
	double G = pCell->GetCellData()->GetItem("G");
	*this->mpOutStream << "G=\"" << target_area << "\" ";
	
	//Perimeter
	 unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    double perimeter = dynamic_cast<VertexBasedCellPopulation<ELEMENT_DIM>*>(pCellPopulation)->rGetMesh().GetSurfaceAreaOfElement(elem_index);
    *this->mpOutStream << "perimeter=\"" << perimeter << "\" ";


	//Neighbours
    std::set<unsigned> neighbours = pCellPopulation->GetNeighbouringLocationIndices(pCell);
	//Number of Neighbours
	int num_neighbours = neighbours.size();
	*this->mpOutStream << "num_neighbours=\"" << num_neighbours << "\" ";
	//List of Neighbouring CellIDs
	*this->mpOutStream << "neighbors=\"";
    for (std::set<unsigned>::iterator neighbour_iter = neighbours.begin();
            neighbour_iter != neighbours.end();
            ++neighbour_iter)
    {
        *this->mpOutStream << " " << (pCellPopulation->GetCellUsingLocationIndex(*neighbour_iter))->GetCellId();
    }
    *this->mpOutStream << "\" ";
	
	//number of Edges	
	VertexElement < ELEMENT_DIM, ELEMENT_DIM > *VertexElement = dynamic_cast<VertexBasedCellPopulation<ELEMENT_DIM>*>(pCellPopulation)->GetElementCorrespondingToCell(pCell);
	int num_edges = VertexElement->GetNumNodes();
	*this->mpOutStream << "num_edges=\"" << num_edges << "\" ";
		
	//End tag   
    *this->mpOutStream << "/>\n";
}

// Dummy implementation because it is required
    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double XMLCellWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return 0.0;
}


// Explicit instantiation
template class XMLCellWriter<1,1>;
template class XMLCellWriter<1,2>;
template class XMLCellWriter<2,2>;
template class XMLCellWriter<1,3>;
template class XMLCellWriter<2,3>;
template class XMLCellWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(XMLCellWriter)
