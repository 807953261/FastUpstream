#ifndef RAPID_JSON_PARSER_H
#define RAPID_JSON_PARSER_H

#include "GiscupParser.h"
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cstring>
#include <cassert>
#include <iostream>
#include "utils.h"


#include "ImportantIdentifier.h"

//Implements our parser that uses RapidJson to parse the JSON file.
template <class GlobalIdType>
class RapidJsonParser: public GiscupParser<GlobalIdType>  {
public:
	typedef  ImportantIdentifier<GlobalIdType> ImportantIdentifierType; //Type of the tokens read from the input JSON...

	//creates the parser. The parser will read the input and create the giscup graph (stored in the base class) 
	RapidJsonParser(const std::string &jsonPath, const std::string &startingPointsPath, bool freeTemporaryMemory = true); //construct parser and parse input.

	void printForDebugging() const;
	~RapidJsonParser() {};

private:	


	//Given the identifiers (having global ids represented as a sequence of characters), create unique
	//"internal" ids for them. I.e., instead of having ids {ABC-123},{EFG-124}, {FFE-999-DC21}, ... we will have 0,1,2,...
	//These "internal" ids are employed to created our graph
	void createInternalIdsForIdentifiers(std::vector<ImportantIdentifierType> &startingPointsList,
                                                       std::vector<ImportantIdentifierType> &controllersList,
                                                        std::vector<array<ImportantIdentifierType,3> > &jsonEdgesList, 
                                                            size_t &numVertices, size_t &numEdges) const;

	//Given the edges, starting points and controllers, create the graph rerpresenting the GISCUP problem
	void createGiscupGraph(const vector< ImportantIdentifierType > &startingPointsIdentifiers, 
                                                const std::vector<array<ImportantIdentifierType,3> > &jsonEdgesList,
                                                const std::vector<ImportantIdentifierType > &jsonControllersList,
                                                int numVertices, int numEdges);

	//Parses the starting point file
	//this implementation is less strict (but slower) than the one available in the FastParallelGiscupParser class (for example, it doesn't require the ids to be surrounded by {})
	void parseStartingPointsFile(std::vector<ImportantIdentifierType > &startingPointsIdentifiers,const std::string &startingPointsPath);
};


#include "RapidJsonParser.cpp"

#endif