/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "RacRhoWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VertexElement.hpp"
#include "VertexBasedCellPopulation.hpp"
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::RacRhoWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("racrho_data.csv")
{
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	
	//AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
	*this->mpOutStream << "# TimeStamp,(Rac, Rho)\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << " ";
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    ///\todo #2645 - if efficiency is an issue, check if this is really needed
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RacRhoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    
	
	
	pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double * rac;
	double * rho;
	rac = new double[num_cells];
	rho = new double[num_cells];
	int i = 0;
	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		double F = cell_iter->GetCellData()->GetItem("F");
		double W = cell_iter->GetCellData()->GetItem("W");
		rac[i] = F;
		rho[i] = W;
		i++;
	}
	for (int j = 0; j<num_cells; j ++){
    	*this->mpOutStream << ",Rac:" << rac[j];
		*this->mpOutStream << ",Rho:" << rho[j];
	}	

}

// Explicit instantiation
template class RacRhoWriter<1,1>;
template class RacRhoWriter<1,2>;
template class RacRhoWriter<2,2>;
template class RacRhoWriter<1,3>;
template class RacRhoWriter<2,3>;
template class RacRhoWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RacRhoWriter)
