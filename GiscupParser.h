#ifndef GISCUP_PARSER_H
#define GISCUP_PARSER_H

#include <cassert>
#include <vector>
#include <string>
#include <map>

using namespace std;

#include "GiscupGraph.h"
#include "ImportantIdentifier.h"

//This class basically represents the GISCUP graph
//It stores the information about the global ids (read from the JSON) of the vertices/edges
//It also stores the graph representing the GICUP problem
//We extend it by adding methods to parse the graph efficiently
template <class GlobalIdType>
class GiscupParser { 
public:
	typedef  ImportantIdentifier<GlobalIdType> ImportantIdentifierType; //Type of the tokens read from the input JSON...

	//given the id of an edge, returns its identifier (we can use this, for example, to output the global id of an edge (given the id we use internally))
	const ImportantIdentifierType& getIdentifierEdge(const int internalId) const {return edgeIdToImportantIdentifier[internalId];} 
	//Same, but for vertices...
	const ImportantIdentifierType& getIdentifierVertex(const int internalId) const {return vertexIdToImportantIdentifier[internalId];}

	const bool isArtificialVertex(const int internalId) const { return internalId >= posStartArtificialVertices; }

	const GiscupGraph &getParsedGraph() const {return graph;}
protected:
	int posStartArtificialVertices; //position of vector vertexIdToImportantIdentifier where the artificial vertices start (all vertices after (including) this position are artificial (representing edges starting points...) )
	std::vector<ImportantIdentifierType > vertexIdToImportantIdentifier;
	std::vector<ImportantIdentifierType > edgeIdToImportantIdentifier;
	bool freeTemporaryMemory;

	GiscupGraph graph;
};

#endif