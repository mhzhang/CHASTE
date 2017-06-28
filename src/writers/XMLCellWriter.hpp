
/*

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/


#ifndef ABSTRACTXMLCELLBASEDWRITER_HPP_
#define ABSTRACTXMLCELLBASEDWRITER_HPP_

#include <string>
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class XMLCellWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
    private:
        std::string mSimulationType;
        std::string mTypeList;
        std::string mAxisDivision;
        std::string mCellCycleModel;
        std::string mExtraSimInfo;

        /** Needed for serialization. */
        friend class boost::serialization::access;
        /**
         * Serialize the object and its member variables.
         *
         * @param archive the archive
         * @param version the current version of this class
         */
        template<class Archive>
            void serialize(Archive & archive, const unsigned int version)
            {
                archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
            }

    public:
        XMLCellWriter();

		/*
        void SetSimulationInfo(std::string simulationType, std::string typeList,
                std::string axisDivision, std::string cellCycleModel,
                std::string extraSimInfo);
		*/

        void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

        void WriteTimeStamp();

        void WriteNewline();

        void OpenOutputFile(OutputFileHandler& rOutputFileHandler);

        double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(XMLCellWriter)

#endif
