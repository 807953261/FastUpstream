#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
#include <cstdio>
#include <sstream>
using namespace std;

//we can enable/disable logging/timming by editing this file
#include "logging.h"
#include "utils.h"
#include "RapidJsonParser.h"
#include "GiscupSolver.h"

#include <set>
using std::set;




#include "GiscupWriter.cpp"



//Solve the GISCUP problem using the RapidJson parser...
//The template parameter defines how we will store the globalIds of the features (we implemented three different ways of doing this)
template <class GlobalIdType>
void solveRapidJsonParser(const string &jsonPath, const string &startingPointsPath, const string &outputPath, bool saveMemory=false) {

	LOGGING(Timer t1("entire main");)
	RapidJsonParser<GlobalIdType> * parser;

	{
		{ //Parse the input, creating a graph
			LOGGING(Timer t("total parsing/reading ");)					
			parser =  new RapidJsonParser<GlobalIdType>(jsonPath,startingPointsPath,saveMemory);
		}
		{
			LOGGING(Timer t("total solve and write results");)
			ofstream out(outputPath.c_str());
			if(!out) {
				clog << "Error creating output file: " << outputPath << endl;
				exit(1);
			}
			vector<int> *outputVerticesIds = new vector<int>();

			//process the graph finding the vertices in simple paths connecting controllers and starting points
			findUpstreamVertices(parser->getParsedGraph(),*outputVerticesIds);

			//given the vertices that should be in the output, find the edges that should be written and write everything to the stream
			writeGlobalIdsOfVerticesAndEdgesInPaths<GlobalIdType, RapidJsonParser<GlobalIdType> >(*parser,parser->getParsedGraph(),*outputVerticesIds, out);
			out.close();
		}
		
		if(saveMemory) {
				delete parser;
		}		
	}
}






//How will the global ids be stored?
const int ID_TYPE_CHAR_PTRS = 0;
const int ID_TYPE_STRINGS = 1;
const int ID_TYPE_GUID = 2;

void solve(const string &jsonPath, const string &startingPointsPath, const string &outputPath, int idType , bool sequential, bool freeMemory) {
	if(sequential)
		omp_set_num_threads(1);

	
	if(idType==ID_TYPE_GUID) { //with rapidJson we can only use two kinds of ways to store the Global IDs... (GlobalIdTypeCharPtrs is not supported)
		solveRapidJsonParser<GlobalIdTypeGUID >(jsonPath, startingPointsPath, outputPath,freeMemory);
	} else {
		solveRapidJsonParser<GlobalIdTypeString >(jsonPath, startingPointsPath, outputPath,freeMemory);
	}
	
}


/*
- The features that are in a simple path between a controller and a starting point are found by computing the block-cut tree
of the graph and, then, iteratively removing leaves without controllers/starting-points (some special cases have to be handled)
- For performance, we implemented a iterative version of Tarjan's algorithm for finding the articulation points of the graph 
(these points are employed for creating the block-cut tree).
- We also tried to parallelize some steps of the algorithm.
- Furthermore, after some benchmarks, we observed that depending on the size of the input files some parsing strategies are more 
suitable than others. Thus, this main.cpp file focus on using an heuristic for selecting a configuration that is fast for the input. 


One particularly important special case happens when a starting point is actually a line:
- If edge e=(a,b) is an starting point --> we create an "artificial/dummy" vertex v and pretend v is an starting point.
- Then, e=(a,b) is replaced with e'=(a,v) and e''=(v,b)   (if a starting point line is composed of several edges --> we do that for all edges)
- This is done during the creation of the graph (in the parsing step) --> thus, our "solver" does not have to handle this special case (the solver
does not have to know that!)

Representing/parsing the data:
- Our parser supports two strategies for storing the 
global ids of the features (from/to/via features describing the connections of the graph):
-- GlobalIdTypeString : the global ids are represented by strings (e.g.: "{3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}"). This has some performance
overheads because we have to store a string with several bytes (for each global id) and, also, because this leads to a lot of
small allocations in the heap.
-- GlobalIdTypeGUID: assumming the global ids always have the format shown above ("{3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}"), we can
codify it using two 64-bit integers (each number/letter is an hexadecimal digit). 

- The available parser was implemented using the RapidJson library (according to the GISCUP website, it is ok to use libraries for assinting
with reading the input files).

The parser throws an integer exception if anything goes wrong --> thus, we start with a faster version and switch to a safer one if something
goes wrong. For example, if we use a parser representing the global ids with the GlobalIdTypeGUID format and during the construction
of an id we detect that it is not in the format we were expecting ("{3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}") --> an exception is thrown
and we switch to another representation. We expect that the input files will be similar to the specification shown in the GISCUP website
but we still do this treatment in case they aren't.
*/

int main(int argc, char **argv) {
	std::ios_base::sync_with_stdio(false); //this makes I/O faster....

	if(argc<4) {
		cerr << "Error: use ./giscup.out inputJsonFilePath startingPointsFilePath outputFilePath" << endl;
		exit(1);
	}

	const string jsonPath = argv[1];
	const string startingPointsPath = argv[2];
	const string outputPath = argv[3];

	//based on the size of the json file we will choose which version of the algorithm to use, the number of threads, etc.
	FILE *fp = fopen(jsonPath.c_str(), "r");
	if(!fp) {
		cerr << "Error: json file cannot be open (is the path correct?): " << jsonPath << endl;
		exit(1);
	}
	size_t jsonFileSizeMB = getFileSize(fp)/(1024*1024);
	fclose(fp);



	LOGGING(clog << "Input json file size: " << jsonFileSizeMB << endl;)

	bool sequential = (jsonFileSizeMB<20); //should we run sequentially?
	bool freeMemory = false;

	//In this contest we tried to squeeze every second..
	//we can save some seconds by not freeing temp. memory (this memory will be fred anyway after the end of the program)
	//also, most of the memory would be needed during the whole lifetime of the program anyway (the diference in the peak
	//memory usage is small)
	freeMemory = jsonFileSizeMB>16000; 

	//calls the solver based on the configuration determined above.
	//if the solver fails (for example, because the global ids are not in a format similar to the expected one) -->
	//we try to use a more generic (safer) solver (that does not make this assumption about the input format)
	
	//In the GISCUP contest we employed a custom parser (that makes some assumptions about the input file)
	//Since the custom parser is not guaranteed to work with input files not holding some assumptions, we decided
	//to not make it available in this version of our program (thus, only the JSON one is available)
	try {
		solve(jsonPath, startingPointsPath, outputPath, ID_TYPE_GUID, sequential, freeMemory);
	} catch(...) {	
		solve(jsonPath, startingPointsPath, outputPath, ID_TYPE_STRINGS, sequential, freeMemory);
	}			
}

