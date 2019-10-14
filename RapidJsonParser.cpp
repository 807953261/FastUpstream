
#define RAPIDJSON_SIMD
#define RAPIDJSON_SSE2


#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
#include <omp.h>
#include <unordered_map>
#include <cstdio>
#include "utils.h"

#include "logging.h"
#include "rapidjson/reader.h"
#include "rapidjson/error/en.h"
#include "rapidjson/filereadstream.h"
#include <condition_variable>
#include <iostream>

using namespace rapidjson;
using namespace std;

//given the start of a global id and the number of characters in it, computes a hash value for it...
//in the FastParallelGiscupParser class this is computed while the input is read
unsigned long long computeHash(const char *start, size_t numChars) {
    unsigned long long hash = 5381;
    int c;                  
    for(size_t i= 0 ;i<numChars;i++,start++) {
        c = *start;
        hash = ((hash << 5) + hash) + c;
    }
    return hash;
}


//This handler deals with the tokens foundin the JSON file
//After reading the input we will have the edges and controllers from the file
template <class ImportantIdentifierType>
struct GiscupJsonHandler : public BaseReaderHandler<UTF8<>, GiscupJsonHandler<ImportantIdentifierType> > {
    GiscupJsonHandler() : insideRows(false), insideControllers(false),foundRowsKeyword(false),foundControllersKeyword(false),next(0),objectImportantInfoRead(0),objectLevel(0), arrayLevel(0) {};
    bool insideRows;
    bool insideControllers;
    bool foundRowsKeyword,foundControllersKeyword;

    //we employ the variables below to keep track about where we are in the json file...
    int objectLevel;
    int arrayLevel;    
    int objectImportantInfoRead;

    //this is the important result of this handler: the edges and controllers from the JSON file..
    vector<array<ImportantIdentifierType,3> > edges;
    vector<ImportantIdentifierType> controllers;

    //temporary tokes employed to store the last "from" token seen during the parsing, the last "to", etc...
    //once we have three of these tokens we can create a new edge!
    ImportantIdentifierType lastFrom, lastTo,lastVia;


    int next; //what is the next kind of token that we will have?
    const static int NEXT_IS_CONTROLLER = 1;
    const static int NEXT_IS_FROM = 2;
    const static int NEXT_IS_TO = 3;
    const static int NEXT_IS_VIA = 4;
    bool String(const char* str, SizeType length, bool copy) {
        switch(next) {
            case NEXT_IS_CONTROLLER:
                controllers.emplace_back(ImportantIdentifierType(-1,str,str+length-1));
                break;
            case NEXT_IS_FROM:  //from
                objectImportantInfoRead++;
                lastFrom = ImportantIdentifierType(-1,str,str+length-1);
                break;
            case NEXT_IS_TO: //to
                objectImportantInfoRead++;
                lastTo = ImportantIdentifierType(-1,str,str+length-1);
                break;
            case NEXT_IS_VIA: //via
                objectImportantInfoRead++;
                lastVia = ImportantIdentifierType(-1,str,str+length-1);
                break;
        }
        if(objectImportantInfoRead==3) {
            edges.push_back({lastFrom,lastTo,lastVia});
            objectImportantInfoRead = 0;
        }

        return true;
    }
    
    bool Key(const char* str, SizeType length, bool copy) { 
        next = 0;
        if (insideRows && objectLevel==2) {
            if(length==12 && strncmp(str,"fromGlobalId",length)==0) { next = NEXT_IS_FROM; } //after this key we will have the global id of a "from" token...
            else if(length==10 && strncmp(str,"toGlobalId",length)==0) { next = NEXT_IS_TO;}
            else if(length==11 && strncmp(str,"viaGlobalId",length)==0) { next = NEXT_IS_VIA;}
        } else if (insideControllers && objectLevel==2) {
            if(length==8 && strncmp(str,"globalId",length)==0) {
                next = NEXT_IS_CONTROLLER;
            }
        } else if(objectLevel==1 && length==4 && strncmp(str,"rows",length)==0) { //we found a "rows" array!
            insideRows = true;
            insideControllers = false;
            objectImportantInfoRead  =0 ;
        } else if(objectLevel==1 && length==11 && strncmp(str,"controllers",length)==0) { //we found a "controllers" array!
            insideRows = false;
            insideControllers = true;
            objectImportantInfoRead  =0 ;
        }
        return true;
    }
    //the edges and controllers info are at the level 2 of the JSON..
    bool StartObject() {  objectLevel++; if(objectLevel==2) { objectImportantInfoRead = 0; } return true; }
    bool EndObject(SizeType memberCount) { objectLevel--; return true; }

    bool StartArray() { return true; }
    bool EndArray(SizeType elementCount) { return true; }
    bool Null() {  return true; }
    bool Bool(bool b) {  return true; }
    bool Int(int i) {  return true; }
    bool Uint(unsigned u) {  return true; }
    bool Int64(int64_t i) { return true; }
    bool Uint64(uint64_t u) {  return true; }
    bool Double(double d) {  return true; }
    bool RawNumber(const char* str, SizeType length, bool copy) {  return true;}
};



//Given the identifiers (having global ids represented as a sequence of characters), create unique
//"internal" ids for them. I.e., instead of having ids {ABC-123},{EFG-124}, {FFE-999-DC21}, ... we will have 0,1,2,...
//These "internal" ids are employed to created our graph
template <class GlobalIdType>
void RapidJsonParser<GlobalIdType>::createInternalIdsForIdentifiers(std::vector<ImportantIdentifierType> &startingPointsList,
                                                       std::vector<ImportantIdentifierType> &controllersList,
                                                        std::vector<array<ImportantIdentifierType,3> > &jsonEdgesList, 
                                                         size_t &numVertices, size_t &numEdges) const {
    LOGGING(Timer t("Creating internal ids for identifiers");)

        //For creating internal ids for the vertices/edges...
    const int numParallelMaps = omp_get_max_threads();

    //this is a "magical" number we employ for tunning the hash table...
    const size_t sizeVertexHashTable = (3*jsonEdgesList.size())/numParallelMaps + 10;
    //we typically have more edges than vertices --> bigger hash table
    const size_t sizeEdgeHashTable   = 2*sizeVertexHashTable;


    //these hash tables map identifiers to internal ids..
    //For example, if we have vertices (v1,v9,v5,v1) --> there are 3 unique vertices and they will be mapped to internal ids 0,1,2
    //For performance, we process the identifiers in parallel using one hash table per thread
    //After filling the tables, we update the vertices ids so that they will all be consistent (and numbered from 0 to numUniqueIdentifiers-1)
    //We label both vertices and edges (but their ids are independent -- we will have a vertex with id 0 and an edge with id 0...)
    unordered_map<ImportantIdentifierType, size_t,typename ImportantIdentifierType::ImportantIdentifierHasherType,typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType> *vertexIdToInternalIdPtr, *edgeIdToInternalIdPtr;
    vertexIdToInternalIdPtr = new unordered_map<ImportantIdentifierType, size_t,typename ImportantIdentifierType::ImportantIdentifierHasherType,typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType> [numParallelMaps];
    edgeIdToInternalIdPtr  = new unordered_map<ImportantIdentifierType, size_t,typename ImportantIdentifierType::ImportantIdentifierHasherType,typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType> [numParallelMaps];
    for(int i=0;i<numParallelMaps;i++) {
        typename ImportantIdentifierType::ImportantIdentifierHasherType hasher;
        typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType comparisor;
        vertexIdToInternalIdPtr[i] = unordered_map<ImportantIdentifierType, size_t,typename ImportantIdentifierType::ImportantIdentifierHasherType,typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType>(sizeVertexHashTable,hasher,comparisor );
        edgeIdToInternalIdPtr[i] = unordered_map<ImportantIdentifierType, size_t,typename ImportantIdentifierType::ImportantIdentifierHasherType,typename ImportantIdentifierType::ImportantIdentifierEqComparisonStringType>(sizeEdgeHashTable,hasher,comparisor );
    }



    {
        LOGGING(Timer t0("Creating ids for the vertices");)
        #pragma omp parallel for schedule(static,1)
        for(int threadId = 0;threadId<numParallelMaps;threadId++) 
            for(auto &edge : jsonEdgesList) {
                {
                    auto &id = edge[0]; //from vertex
                    const int bucketToUse = id.hash%numParallelMaps;
                    if(bucketToUse == threadId) { //siince each thread uses a different bucket, we may insert into the hash tables in parallel
                        const auto internalIdIt = vertexIdToInternalIdPtr[bucketToUse].find(id);
                        if(internalIdIt ==vertexIdToInternalIdPtr[bucketToUse].end()) {
                            id.id = vertexIdToInternalIdPtr[bucketToUse].size();
                            vertexIdToInternalIdPtr[bucketToUse][id] = id.id;
                        } else {
                            id.id = internalIdIt->second;
                        }
                    }
                }
                {
                    auto &id = edge[1]; //to vertex
                    const int bucketToUse = id.hash%numParallelMaps;
                    if(bucketToUse == threadId) {
                        const auto internalIdIt = vertexIdToInternalIdPtr[bucketToUse].find(id);
                        if(internalIdIt ==vertexIdToInternalIdPtr[bucketToUse].end()) {
                            id.id = vertexIdToInternalIdPtr[bucketToUse].size();
                            vertexIdToInternalIdPtr[bucketToUse][id] = id.id;
                        } else {
                            id.id = internalIdIt->second;
                        }
                    }
                }
            }        
    }

    {
        //set the ids of the global ids (controllers...)
        #pragma omp parallel for schedule(static,1)
        for(int threadId = 0;threadId<numParallelMaps;threadId++) {
            for(ImportantIdentifierType &id : controllersList) { 
                const int bucketToUse = id.hash%numParallelMaps;
                if(bucketToUse == threadId) {
                    const auto internalIdIt = vertexIdToInternalIdPtr[bucketToUse].find(id);
                    if(internalIdIt ==vertexIdToInternalIdPtr[bucketToUse].end()) {
                        id.id = vertexIdToInternalIdPtr[bucketToUse].size();
                        vertexIdToInternalIdPtr[bucketToUse][id] = id.id;
                    } else {
                        id.id = internalIdIt->second;
                    }
                }
            }
        }

    }

    //Now each hash table will have the ids of the vertices in them... however, different vertices from different hash tables may have the same ids
    //Ex: vertexIdToInternalIdPtr[0] will have a vertex with id 0, vertexIdToInternalIdPtr[1] will also have a (different) vertex with id 0...
    //We fix this by adding the size of the previous hash tables to the ids of the vertices.
    //Ex: if vertexIdToInternalIdPtr[0] has 100 vertices, their vertices ids will be 0,1,...99 and if vertexIdToInternalIdPtr[1] has 50 elements,
    //the ids will be 100 + 0, 100 + 1, 100 + 2, ... 100+49 (100 is the sum of the sizes of the previous hash tables)
    vector<size_t> startIdsForVerticesEachBucket(numParallelMaps,0);
    { 
        LOGGING(Timer t0("aggregating the vertice ids from the different hash tables");)
        for(int i=1;i<numParallelMaps;i++) 
            startIdsForVerticesEachBucket[i] =  startIdsForVerticesEachBucket[i-1] + vertexIdToInternalIdPtr[i-1].size();

        const size_t numEdges = jsonEdgesList.size();
        #pragma omp parallel for
        for(size_t i=0;i<numEdges;i++) {
            {
                ImportantIdentifierType &id = jsonEdgesList[i][0]; //from
                int bucketToUse = id.hash%numParallelMaps;
                id.id =  id.id + startIdsForVerticesEachBucket[bucketToUse];                
            }
            {
                ImportantIdentifierType &id = jsonEdgesList[i][1]; //to
                int bucketToUse = id.hash%numParallelMaps;
                id.id =  id.id + startIdsForVerticesEachBucket[bucketToUse];                
            }
        }               
        
        //update the ids of the global ids (controllers...)
        const size_t numIdsThisVector = controllersList.size();
        #pragma omp parallel for
        for(size_t i=0;i<numIdsThisVector;i++) {
            ImportantIdentifierType &id = controllersList[i];
            const int bucketToUse = id.hash%numParallelMaps;
            id.id =  id.id + startIdsForVerticesEachBucket[bucketToUse];                  
        }               
        
    }

    //Creating the via ids...
    {
        LOGGING(Timer t0("creating the internal ids for the via tokens");)
        #pragma omp parallel for schedule(static,1)
        for(int k=0;k<numParallelMaps;k++) 
            for (auto &edge:jsonEdgesList)  {
                ImportantIdentifierType &id = edge[2];
                int bucketToUse = id.hash%numParallelMaps;
                if(bucketToUse!=k) continue;
                
                const auto internalIdIt = edgeIdToInternalIdPtr[bucketToUse].find(id);
                if(internalIdIt ==edgeIdToInternalIdPtr[bucketToUse].end()) {
                    id.id = edgeIdToInternalIdPtr[bucketToUse].size();
                    edgeIdToInternalIdPtr[bucketToUse][id] = id.id;
                } else {
                    id.id = internalIdIt->second;
                }                   
            }
                       
    }
    vector<size_t> startIdsForEdgesEachBucket(numParallelMaps,0);
    { 
        LOGGING(Timer t0("aggregating the edge (via) ids from the different hash tables");)
        for(int i=1;i<numParallelMaps;i++) 
            startIdsForEdgesEachBucket[i] =  startIdsForEdgesEachBucket[i-1] + edgeIdToInternalIdPtr[i-1].size();

        size_t numEdges = jsonEdgesList.size();
        #pragma omp parallel for
        for(size_t j=0;j<numEdges;j++) {        
            auto &edge = jsonEdgesList[j];
                        
            ImportantIdentifierType &id = edge[2];
            
            int bucketToUse = id.hash%numParallelMaps;
            id.id =  id.id + startIdsForEdgesEachBucket[bucketToUse];                       
        }
    }



    //Creates ids for the starting points...
    size_t numStartingPointIdentifiers = startingPointsList.size();
    #pragma omp parallel for
    for(size_t j=0;j<numStartingPointIdentifiers;j++) {
        ImportantIdentifierType &i = startingPointsList[j];
        int bucketToUse = i.hash%numParallelMaps;
        const auto vertexInternalIdIt = vertexIdToInternalIdPtr[bucketToUse].find(i);
        if(vertexInternalIdIt!=vertexIdToInternalIdPtr[bucketToUse].end()) { //is it a vertex?
            i.id = vertexInternalIdIt->second + startIdsForVerticesEachBucket[bucketToUse];
            i.setType(ImportantIdentifierType::VERTEX_STARTING_POINT_IDENTIFIER);
        } else { //if it is not an edge --> it doesn't exist in the graph!
            const auto edgeInternalIdIt = edgeIdToInternalIdPtr[bucketToUse].find(i);
            if(edgeInternalIdIt==edgeIdToInternalIdPtr[bucketToUse].end()) {                
                i.id=-1;
                continue;
            } 
            i.id = edgeInternalIdIt->second  + startIdsForEdgesEachBucket[bucketToUse];
            i.setType(ImportantIdentifierType::EDGE_STARTING_POINT_IDENTIFIER);
        }
    }



    numVertices = 0; // vertexIdToInternalId.size();
    numEdges = 0;
    for(int i=0;i<numParallelMaps;i++) {
        numVertices += vertexIdToInternalIdPtr[i].size();
        numEdges += edgeIdToInternalIdPtr[i].size();
    } 


    if(GiscupParser<GlobalIdType>::freeTemporaryMemory) {
        LOGGING(Timer t("Freeing parser temporary memory");)
        delete []vertexIdToInternalIdPtr;
        delete []edgeIdToInternalIdPtr;
    }


}





//given the identifiers with internal integer ids (i.e., ids from 0...nv-1 for vertices, 0...ne-1 for edges),
//create the graph (adjacency list) representing the problem
template <class GlobalIdType>
void RapidJsonParser<GlobalIdType>::createGiscupGraph(const vector< ImportantIdentifierType > &startingPointsIdentifiers, 
                                                const std::vector<array<ImportantIdentifierType,3> > &jsonEdgesList,
                                                const std::vector<ImportantIdentifierType > &jsonControllersList,
                                                int numVertices, int numEdges) {

    LOGGING(Timer t("Creating graph and auxiliary information");)

    LOGGING(cerr << "Creating graph with: " << numVertices << " vertices" << endl;)    
    

    LOGGING(cerr << "Creating edges... with " << numEdges << endl;)
    GiscupParser<GlobalIdType>::edgeIdToImportantIdentifier.resize(numEdges);

    LOGGING(cerr << "after edge id to important ids... " << endl;)

    //isStartingPointE[i] = true iff the edge with internal id i is s starting point...
    vector<int> isStartingPointE(numEdges,false);



    //label the starting points Edges starting points...
    //maybe do this in parallel
    #pragma omp parallel
    for (const ImportantIdentifierType &i:startingPointsIdentifiers) {
        if(i.id==-1) continue; //maybe this vertex does not exist in the graph (but exists in the list of starting points...)
        if(i.type==ImportantIdentifierType::EDGE_STARTING_POINT_IDENTIFIER){
            //edgeStartingPoints.push_back(i.id);
            isStartingPointE[i.id] = true;
        }
    }


    size_t numArtificialVertices = 0;

    //compute the degree of the regular (non artificial) vertices
    //also counts the number of edges that are starting point --> these edges will lead to the creation of artificial vertices
    vector<size_t> regularVertexDegrees(numVertices,0);
    {
        LOGGING(Timer t2("count number of artificial vertices");)

        //#pragma omp parallel for reduction(+:numArtificialVertices)
        for(const auto &edge:jsonEdgesList){
            regularVertexDegrees[edge[0].id]++; //from
            regularVertexDegrees[edge[1].id]++; //to  
            if(isStartingPointE[edge[2].id]) {
                numArtificialVertices++;
            }           
        }        
    }
   

    //we will have to create at least numStartingPointEdges artificial vertices to represent edge starting points...
    this->posStartArtificialVertices =numVertices;
    
    GiscupParser<GlobalIdType>::vertexIdToImportantIdentifier.resize(numVertices+numArtificialVertices); //given the internal id of an identifier, gets the identifier.
   
    GiscupParser<GlobalIdType>::graph  = std::move(GiscupGraph(numVertices+numArtificialVertices,regularVertexDegrees,2) ); //adjacency list...

    //label the starting points as starting points...
    #pragma omp parallel
    for (const ImportantIdentifierType &i:startingPointsIdentifiers) {
        if(i.id==-1) continue; //maybe this vertex does not exist in the graph (but exists in the list of starting points...)
        if(i.type==ImportantIdentifierType::VERTEX_STARTING_POINT_IDENTIFIER) {
            GiscupParser<GlobalIdType>::graph.setStartingPointVertex(i.id,true);
        } 
    }
    LOGGING(cerr << "Creating adjacency information" <<endl;)
    

    size_t ctStartingPointEdgesSoFar = 0;
    {
        LOGGING(Timer t2("Setting controllers info...");)

        for(const ImportantIdentifierType &identifier:jsonControllersList)  {
            if(identifier.id==-1) continue; //ignore controllers that do not exist in the graph!
            GiscupParser<GlobalIdType>::graph.setControllerVertex(identifier.id,true); 
        }       
    }

    LOGGING(Timer t2("Second part of graph creation...");)
    //create the adjacency list and the list of controllers



    //Maybe a controller
    #pragma omp parallel for
    for(size_t i =0;i<jsonControllersList.size();i++) {
        GiscupParser<GlobalIdType>::vertexIdToImportantIdentifier[jsonControllersList[i].id] = jsonControllersList[i];
    }

    for(const auto &edge: jsonEdgesList) {        
        //if the "via" is a vertex --> we ignore this connectivity!
        //it only happens when a point connect to itself via itself 
        //this does not change anything in the results...
        const ImportantIdentifierType &from = edge[0];
        const ImportantIdentifierType &to = edge[1];
        const ImportantIdentifierType &via = edge[2];

        GiscupParser<GlobalIdType>::vertexIdToImportantIdentifier[from.id] = from;
        GiscupParser<GlobalIdType>::vertexIdToImportantIdentifier[to.id] = to;

        
        GiscupParser<GlobalIdType>::edgeIdToImportantIdentifier[via.id] = via; //can the via be a vertex?
    
        if(isStartingPointE[via.id]) {  
            //if the edge is a starting point --> we now create an artificial vertex...

            ImportantIdentifierType &viaArtificialVertex = GiscupParser<GlobalIdType>::vertexIdToImportantIdentifier[this->posStartArtificialVertices+ctStartingPointEdgesSoFar];             

            viaArtificialVertex = (GiscupParser<GlobalIdType>::getIdentifierEdge(via.id));

            viaArtificialVertex.id = this->posStartArtificialVertices+ctStartingPointEdgesSoFar;
            ctStartingPointEdgesSoFar++;
            
            //viaArtificialVertex.id = vertexIdToImportantIdentifier.size()-1;
            GiscupParser<GlobalIdType>::graph.setStartingPointVertex(viaArtificialVertex.id,true); //isStartingPointV[viaArtificialVertex.id] = true;
            GiscupParser<GlobalIdType>::graph.setControllerVertex(viaArtificialVertex.id,false); //isControllerV[viaArtificialVertex.id] = false;             

            
            GiscupParser<GlobalIdType>::graph.addEdge(from.id,viaArtificialVertex.id,via.id);
            GiscupParser<GlobalIdType>::graph.addEdge(viaArtificialVertex.id,from.id,via.id);

            GiscupParser<GlobalIdType>::graph.addEdge(to.id,viaArtificialVertex.id,via.id);
            GiscupParser<GlobalIdType>::graph.addEdge(viaArtificialVertex.id,to.id,via.id);
        } else {
            //loop : we only deal with a loop if it is a starting point (otherwise we ignore them since they will never be in the output...)
            if(from.id!=to.id) {
                GiscupParser<GlobalIdType>::graph.addEdge(from.id,to.id,via.id);
                GiscupParser<GlobalIdType>::graph.addEdge(to.id,from.id,via.id);
            }
        }                       
        
    }
    
}




template <class GlobalIdType>
void RapidJsonParser<GlobalIdType>::parseStartingPointsFile(std::vector<ImportantIdentifierType > &startingPointsIdentifiers,const std::string &startingPointsPath)  {
    LOGGING(Timer t("finding starting point identifiers");)

    ifstream fin(startingPointsPath.c_str());
    string startingPoint;
    while(getline(fin,startingPoint)) {
        //trim the string...
        startingPoint.erase(startingPoint.begin(),std::find_if(startingPoint.begin(), startingPoint.end(),std::not1(std::ptr_fun<int,int>(std::isspace))));
        startingPoint.erase(std::find_if(startingPoint.rbegin(),startingPoint.rend(),std::not1(std::ptr_fun<int, int>(std::isspace))).base(), startingPoint.end());
        if(startingPoint.size()==0) continue; //ignore empty ids...

        char * stringStart = &(startingPoint[0]);
        startingPointsIdentifiers.push_back(ImportantIdentifierType(-1,stringStart,stringStart+startingPoint.size()-1));
    }
}



//Creates the parser and parses the json/starting point files
template <class GlobalIdType>
RapidJsonParser<GlobalIdType>::RapidJsonParser(const std::string &jsonPath, const std::string &startingPointsPath, bool freeTemporaryMemory) {
    this->freeTemporaryMemory = freeTemporaryMemory;
    GiscupJsonHandler<ImportantIdentifier<GlobalIdType> > giscupHandler;
    
    std::vector<ImportantIdentifierType> &controllersList = giscupHandler.controllers; 
    std::vector<array<ImportantIdentifierType,3>> &jsonEdgesList = giscupHandler.edges; 
    { //Use rapidJson to parse the json file. For performance, instead of constructing the entire JSON tree we use a parser based on events...
        LOGGING(Timer t0("Parse json");)   
        Reader reader;
        size_t buffSize = ((size_t)64)*1024*1024;
        char *readBuffer= new char[buffSize];
        FILE* fp = fopen(jsonPath.c_str(), "r"); 
        FileReadStream is(fp, readBuffer, buffSize);        
        reader.Parse<rapidjson::kParseNumbersAsStringsFlag>(is,giscupHandler); //we do not need the numbers -- we can parse than as strings for performance (thus, the parse will not try to convert them to integer/double values)
        fclose(fp);
        if(freeTemporaryMemory) delete []readBuffer;
    }
    std::vector<ImportantIdentifierType> startingPointsList; 
    parseStartingPointsFile(startingPointsList,startingPointsPath);

    //Create the internal ids (zero-based) to represent the global ids
    size_t numVertices = 0;
    size_t numEdges = 0;
    {
        LOGGING(Timer t0("creating internal ids for global ids");)
        createInternalIdsForIdentifiers(startingPointsList,controllersList,jsonEdgesList,numVertices,numEdges);
    }

    //Finally, create the graph representing the problem
    createGiscupGraph(startingPointsList, jsonEdgesList,controllersList,numVertices,numEdges);

}