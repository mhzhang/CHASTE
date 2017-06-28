
/*

Author: MoHan Zhang
Contact: mohan_z@hotmail.com

*/

#ifndef SHAPEWRITER_HPP_
#define SHAPEWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing the cell population
 * adjacency (i.e. connectivity) matrix to file.
 *
 * The output file is called cellpopulationadjacency.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ShapeWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
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
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    ShapeWriter();
	
	
	void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the population to visit.
     */
    void VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
	
	
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ShapeWriter)

#endif /* CELLPOPULATIONADJACENCYMATRIXWRITER_HPP_ */
