#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <map>
#include <cstdio>
#include <stack>
#include <algorithm>
using namespace std;

#include "logging.h"
#include "utils.h"
#include "GiscupSolver.h"





//Implementation of the classical Tarjan's algorithm (for finding articulation points/biconnected components of a graph)
//It works even for disconnected graphs, multigraphs and non-simple graphs
//We could've used Boost's implementation, but this one is faster.
//This algorithm is typically recursive -- but we implemented it using a stack to simulate the recursion (this prevents stack overflows that 
//could happen in big graphs)
//Finds articulation points and biconnected components in a graph


//ATTENTION:
//For performance, since the entire graph is traversed anyway --> we also check if the connected components (CCs) contain both kinds of important vertices
//(controllers and starting points). If a connected component does not have both --> all its biconnected components and articulations
//are ignored...

//Furthermore, since we have to access each vertex at least once --> we also fill the verticesBothControllersStartingPoints vector
//verticesBothControllersStartingPoints contains all vertices of the graph that are both a controller and a starting point
//we could have done this outside this function, but we do it here to save some time...


//Input: a graph represented by an adjacency list (GiscupGraph)
//Output:
//compRaggedArray and compRaggedArrayStart: ragged array (for performance) of biconnected components of the graph 
//each "subarray" in compRaggedArray is a biconnected component. 
//Subarray 0 starts at position compRaggedArrayStart[0] of compRaggedArray
//Subarray 1 starts at position compRaggedArrayStart[1] of compRaggedArray
// ....
//The last position of compRaggedArrayStart represents the position the after last biconnected component would start
//(i.e., the number of biconnected components is compRaggedArrayStart.size()-1)

//articulations: vector of indices of the articulation points of the graph
//verticesBothControllersStartingPoints: vector of indices of the vertices of the graph that are both a controller and a starting point
void biconnectedComponents_ignoringCCsWOImportantVertices(const GiscupGraph &graph,
														vector<int> &compRaggedArray,vector<int> &compRaggedArrayStart,
														vector<int> &isArticulation, vector<int> &articulations, vector<int> &verticesBothControllersStartingPoints) {

	
	LOGGING(Timer t("Run non-recursive tarjan");)
	const int n = graph.numVertices();

	compRaggedArray.reserve(n);
	


	vector<pair<int,int> > numLow(n,pair<int,int>(-1,-1));

	vector< int > vertexStack; 
	isArticulation.clear();
	isArticulation.resize(n,0);

	struct CallStackState {
		CallStackState(int currentChildPos_, int prev_, int parent_): parent(parent_),prev(prev_),currentChildPos(currentChildPos_) {}
		int currentChildPos, prev, parent;
	};

	vector< CallStackState > callStack; //we use this to simulate a call stack (so that we can implement the algorithm non recursivelly).
																		//array[0] is the next child of the current node that we have to process..
																		//array[1] is the previous node (before parent...)
																		//array[2] is the parent node
																		
 	vertexStack.reserve(n);
 	callStack.reserve(n);


	for(int i=0;i<n;i++) {
		//if the vertex is a controller and a starting point --> it will be in the output! (we treat this as a special case of our algorithm)
		//we employ artificial vertices to represent starting points that are edges --> (but an artificial vertex will never be a controller...)
		//we could have treated this somewhere else --> but here we will traverse all the vertices anyway!
		if(graph.isControllerVertex(i) && graph.isStartingPointVertexOrArtVertex(i)) verticesBothControllersStartingPoints.push_back(i);

		if(numLow[i].first!=-1) continue; //have I already visited this node?
		vertexStack.clear();
		vertexStack.push_back(i);
		bool foundControllerInCC = graph.isControllerVertex(i);
		bool foundStartingPointInCC = graph.isStartingPointVertexOrArtVertex(i);
		int oldArticulationsSize = articulations.size();
		int oldCompRaggedArraySize = compRaggedArray.size();
		int oldCompRaggedArrayStartSize = compRaggedArrayStart.size();

		callStack.clear();
		callStack.emplace_back(0,i,i);
		int numChildrenRoot = 0;
		numLow[i].first = numLow[i].second = 0;

		int currentLabel = 1;
		while(!callStack.empty()) {
			int &currentChildPos = callStack.back().currentChildPos;
			int prev = callStack.back().prev;
			int parent = callStack.back().parent;			
			
			const int numNeighbors = graph.vertexDegree(parent);
			if(currentChildPos<numNeighbors) {
				int child = graph.getAdjVertex(parent,currentChildPos);
				currentChildPos++;

				if(prev == child) continue; //edge pointing to previous node..
				if(numLow[child].first==-1) { //I have never visited this node...
					vertexStack.emplace_back(child);
					foundControllerInCC |= graph.isControllerVertex(child);
					foundStartingPointInCC |= graph.isStartingPointVertexOrArtVertex(child);

					callStack.emplace_back(0,parent,child);
					numLow[child].first = numLow[child].second = currentLabel++;					
				} else if(numLow[child].first <= numLow[parent].first){
					if(numLow[child].first<numLow[parent].second) numLow[parent].second = numLow[child].first;
				}
			}
			else {
				callStack.pop_back(); //I have already processed all children of parent...
				if(callStack.size()>1) {
					if(numLow[parent].second >= numLow[prev].first) {
						//the edges from the position containing (prev,curr) is a new biconnected component...
						compRaggedArrayStart.push_back(compRaggedArray.size());
						compRaggedArray.push_back(prev);

						int i = vertexStack.size()-1;
						if(prev!=parent) {
							while(vertexStack[i]!=parent) i--;
							compRaggedArray.insert(compRaggedArray.end(),vertexStack.begin()+i,vertexStack.end());
							vertexStack.resize(i);
						}

						if(!isArticulation[prev]) {
							isArticulation[prev] = true;
							articulations.emplace_back(prev);
						}
					}
					if(numLow[parent].second<numLow[prev].second) numLow[prev].second = numLow[parent].second;
				} else if(!callStack.empty()) {
					numChildrenRoot++;

					compRaggedArrayStart.push_back(compRaggedArray.size());
					compRaggedArray.push_back(prev);

					int i = vertexStack.size()-1;
					if(prev!=parent) {
						while(vertexStack[i]!=parent) i--;
						compRaggedArray.insert(compRaggedArray.end(),vertexStack.begin()+i,vertexStack.end());
						vertexStack.resize(i);
					}
				}
			}
		}
		

		if(numChildrenRoot>1) {
			//root is an articulation point...
			if(!isArticulation[i]) {
				isArticulation[i] = true;
				articulations.push_back(i);
			}
		}

		//ignore connected components (CCs) not containing both kinds of important vertices...
		//this may happen in graphs with multiple connected components (or even in graphs without starting points/controllers)
		if(!foundControllerInCC || !foundStartingPointInCC) {
			//this connected component does not have both a controller and a starting point... we will not add it to the ragged array...
			compRaggedArray.resize(oldCompRaggedArraySize);
			articulations.resize(oldArticulationsSize);
			compRaggedArrayStart.resize(oldCompRaggedArrayStartSize);
		}
	}
	//this is the starting position of the "after the last" biconnected component
	//(we have this "dummy" position because it makes easier to get the size of the components: size[i] = compRaggedArrayStart[i+1]-compRaggedArrayStart[i] )
	//(without having to check if i is the last component)
	compRaggedArrayStart.push_back(compRaggedArray.size()); 
}





//Represents a vertex of the block-cut tree
//Each vertex represents a block (biconnected component) (or an articulation point) of the original graph
//The vertices have some flags/information and the ids of the neighbors 
struct BlockCutTreeVertex {
	BlockCutTreeVertex(): disabled(false),vertexInfo(0) {}
	inline void setFlagHasController() {  (vertexInfo|=1); }
	inline void setFlagHasStartingPointVertex() {  (vertexInfo|=2); }
	inline void setFlagHasMultipleImportantVertices() { (vertexInfo|=4); }
	inline void setFlagDisabled() { disabled = true; }


	inline bool flagDisabled() const  {return disabled;}
	inline bool flagHasImportantVertex() const  { return (vertexInfo!=0); }
	inline bool flagBlockWithoutImportantVertex() const  { return (vertexInfo==0); }
	inline bool flagHasController() const  { return (vertexInfo&1); }
	inline bool flagHasStartingPointVertex() const  { return (vertexInfo&2); }
	inline bool flagHasMultipleImportantVertices() const  { return (vertexInfo&4); }
	 
	inline bool isLeaf() const {return numActiveNeighbors<=1;}

	int numActiveNeighbors; //number of neighbors that have not been removed (disabled) from the tree yet
	short vertexInfo;
	bool disabled; //have we removed this vertex from the BC-tree?
	vector<int> neighbors; //ids of the neighbors of this vertex...
};


//Creates a block-cut tree
//returns the adjacency information of each block (stored in bcTreeVertices)
//This function also works for disconnected graphs (in this case it will create a forest of BC-trees)

//each vertex (stored in bcTreeVertices) representing a block (or articulation) stores some information (filled by this function):
//has it been disabled: we disable vertices when we "remove" it from the tree (for performance, we actually do not remove -- we only set as disabled)
//how many active neighbors the vertex has (initially this is equal to the number of neighbors)
//does it have an important vertex? (articulation or starting point)
//does it have a controller vertex?
//does it have a starting point ?
//does it have multiple important vertices?

//each block represents a biconnected component or an articulation point
//blockId is also filled by this function
//blockId is the block of each graph vertex (articulations are in multiple block --> we assume their block is their own tree block)
//numEdgesTree is the number of edges created in the BC-Tree
void createBCTree(int &numEdgesTree, vector<BlockCutTreeVertex > &bcTreeVertices,vector<int> &blockId,  
					const vector<int> &biconnectedCompsRaggedArray, const vector<int> &biconnectedCompsRaggedArrayStart, const vector<int> &isArticulation, const vector<int> &articulations, const GiscupGraph &graph) {
	LOGGING(Timer t("Create BC tree");)
	//in the BC tree, there is one vertex for each articulation and for each block...
	
	int numBiconnectedComponents = biconnectedCompsRaggedArrayStart.size()-1;
	int nvTree = numBiconnectedComponents + articulations.size();

	blockId.clear();
	blockId.resize(graph.numVertices()); //id of the tree vertex containing each graph vertex....

	bcTreeVertices.clear();
	bcTreeVertices.resize(nvTree);
	numEdgesTree = 0;

	#pragma omp parallel for
	for(unsigned i=0;i<articulations.size();i++) {
		blockId[articulations[i]] = i; //the first articulations.size() blocks will represent articulations...

		if(graph.isControllerVertex(articulations[i])) bcTreeVertices[i].setFlagHasController();
		if(graph.isStartingPointVertexOrArtVertex(articulations[i])) bcTreeVertices[i].setFlagHasStartingPointVertex();
	}


	vector<int> articulationDegreesApprox(articulations.size(),0);
	
	//Now, we create the BC Tree, making each biconnected component neighbor of the articulations in it
	#pragma omp parallel
	{
		vector<int> myArts; //to avoid vector reallocations, we use this vector to tempraly copy the articulations in this block
												//after the for loop, we insert these articulations in the adjacency list (insertion is faster than many push backs...)
		myArts.reserve(100);

		#pragma omp for 
		for(unsigned i=0;i<numBiconnectedComponents;i++) {
			int compId = articulations.size()+i; //id (index in the bcTreeVertices) of this biconnected component in the tree
			int compStart = biconnectedCompsRaggedArrayStart[i];
			int compEnd = biconnectedCompsRaggedArrayStart[i+1];

			myArts.resize(0); //articulations in this biconnected component...
			for(int i=compStart;i<compEnd;i++) { //now we process the vertices (from the graph original) inside this biconnected component
				int v = biconnectedCompsRaggedArray[i]; //id of the vertex in the original graph
				if(!isArticulation[v]) {
					blockId[v] = compId; //the block id of a biconnected component in the tree will be equal to its id + the number of articulations
				} else {
					int articulationId = blockId[v];
					myArts.push_back(articulationId); //in the BC-Tree, each block will have an edge to the block of each of its articulations.
				}
				if(graph.isImportantVertex(v)) {
					//we set "has multiple important vertices" when we find a second important vertex in the tree block
					if(bcTreeVertices[compId].flagHasImportantVertex()) bcTreeVertices[compId].setFlagHasMultipleImportantVertices();
					if(graph.isControllerVertex(v)) bcTreeVertices[compId].setFlagHasController();
					if(graph.isStartingPointVertexOrArtVertex(v)) bcTreeVertices[compId].setFlagHasStartingPointVertex();
				}				
			}

			for(int a:myArts) articulationDegreesApprox[a]++; //this is not thread-safe, but it is ok (we only need an approximation to reserve memory and accelerate the vector push-backs) - we observed this is faster than using atomics!
			bcTreeVertices[compId].neighbors.reserve(myArts.size());
			//Each biconnected component in the tree is neighbor of the articulations in it...
		  //later we will make the articulations neihbor of the biconnected components
			bcTreeVertices[compId].neighbors.insert(bcTreeVertices[compId].neighbors.end(),myArts.begin(),myArts.end());
		}		
	}

	//Now, add an edge (to the BC-tree) from the articulations to the BC-tree vertices representing the biconnected components having the articulation...
	for(int i=0;i<articulations.size();i++)
			bcTreeVertices[i].neighbors.reserve(articulationDegreesApprox[i]); //this is just an approximation for the size...

	for(unsigned i=0;i<numBiconnectedComponents;i++) {
		int compId = articulations.size()+i;

		int numEdges=  bcTreeVertices[compId].neighbors.size();
		for(int i=0;i<numEdges;i++) {
			int articulationId = bcTreeVertices[compId].neighbors[i];
			bcTreeVertices[articulationId].neighbors.push_back(compId); //if the articulation is neighbor of a Biconnected Component --> the bic. comp. is neighbor of hte articulation
		}
		numEdgesTree+= numEdges;
	}
	
}



//Given a graph (represented by a GiscupParser), computes the list of vertices that are in a path between a controller and a starting point
//The graph must have the starting point and controllers vertices labeled
//This function does not treat the special case when the starting point is an edge
//HOWEVER--> we handle this by creating artificial vertices representing these starting point edges (this is handled before this function is called)
//this greatly simplify the treatment of the special cases...
void findUpstreamVertices(const GiscupGraph &graph,vector<int> &outputVerticesIds) {
	LOGGING(Timer t0("find upstream vertices function"));
	
	//biconnected components... (represented using a ragged array for performance)
	vector<int> biconnectedCompsRaggedArray,biconnectedCompsRaggedArrayStart;

	vector<int> isArticulation;
	vector<int> articulations;


	int numConnectedComponentsGraph;
	//Use Tarjan's algorithm for finding the biconnected components.
	//Since this algorithm will process all the vertices anyway --> it also fills the "isArticulation" vector
	//and it also adds to the output the vertices that are simultaneously controllers and starting points (could this happen? I think
	//this is a legal input...) 
	biconnectedComponents_ignoringCCsWOImportantVertices(graph,biconnectedCompsRaggedArray,biconnectedCompsRaggedArrayStart,isArticulation,
																			articulations,outputVerticesIds);


	LOGGING(cerr << "Number of biconnected components found: " << biconnectedCompsRaggedArrayStart.size()-1 << "\n";)
	LOGGING(cerr << "Number of articulations found: " << articulations.size() << "\n";)



	vector<BlockCutTreeVertex> bcTreeVertices; //graph representing the Block-cut tree
	vector<int> blockId; //blockId is the block of each vertex (articulations are in multiple block --> we assume their block is their own tree block)


	int numEdgesBcTree;
	//Create the block-cut tree (BC-Tree) of the graph
	createBCTree(numEdgesBcTree,bcTreeVertices,blockId, biconnectedCompsRaggedArray,biconnectedCompsRaggedArrayStart,isArticulation, articulations, graph);
	int nvTree = bcTreeVertices.size();
	LOGGING(cerr << "Number of vertices in BC-tree: " << nvTree << "\n";)
	LOGGING(cerr << "Number of edges in BC-tree: " << numEdgesBcTree << "\n";)



	//process the tree, removing leaves w/o important vertices (leaves without starting-point/controllers cannot be in a simple path between a controller and a starting point)...
	//This process is repeated because a vertex that becomes a new leaf can also be removed!
	//We actually do not "remove" the leaves --> we label them as "disabled" and update the numActiveNeighbors counter (which simulates the degree of each tree vertex...) 
	{
		LOGGING(Timer t("Removing leaves from tree");)
		stack<int> leavesTreeToRemove;		
		for(int v=0;v<nvTree;v++) {
			bcTreeVertices[v].numActiveNeighbors = bcTreeVertices[v].neighbors.size(); 
			if(bcTreeVertices[v].isLeaf() && bcTreeVertices[v].flagBlockWithoutImportantVertex()) {
				leavesTreeToRemove.push(v);
			}
		}
		//Leaves have in-degree (degree, since the graph is not directed) <= 1 (0 if the leaf is a single connected component)

		int numLeavesRemoved = 0;
		while(!leavesTreeToRemove.empty()) {
			int v = leavesTreeToRemove.top(); leavesTreeToRemove.pop();

			numLeavesRemoved++;
			bcTreeVertices[v].setFlagDisabled(); //"remove" the leaf
			
			//updates the degree of the neihbors
			for(int neighbor:bcTreeVertices[v].neighbors) {
				bcTreeVertices[neighbor].numActiveNeighbors--;

				//for each non-disabled neighbor that is a new leaf (and does not contain an important vertex) --> it will be removed...
				if(!bcTreeVertices[neighbor].flagDisabled() && bcTreeVertices[neighbor].isLeaf() && bcTreeVertices[neighbor].flagBlockWithoutImportantVertex())
					leavesTreeToRemove.push(neighbor);
			}
		}
		LOGGING(cerr << "Number of vertices in the BC-Tree and number of leaves removed: " << nvTree << " " << numLeavesRemoved << endl;)
	}




	//Now, we have a BC-Tree (actually a forest -- if the graph is disconnected) where all leaves (unless the tree is empty)
	//have important vertices (starting points or controllers). We also know that all the connected components have both
	//controllers and starting points --> all the vertices in the remaining biconnected components should be in a path between a controller and a starting point

	//SPECIAL CASE: however, if the only important vertex of a leaf is an articulation still in the tree --> the other vertices in this leave
	//will never be in a simple path between a controller and a starting point.

	{
		LOGGING(Timer t("select vertices to write to output");)
		//the vertices of the tree associated to the articulations do not need to be added to the output... (there will be at least a copy of each of them in one of the blocks...)
		
		for(int tv=articulations.size();tv<nvTree;tv++) { //for each tree vertex block (that is not an articulation)
			if(bcTreeVertices[tv].flagDisabled()) continue; //ignore removed blocks...
			int idComponent = tv - articulations.size(); //id of the biconnected component associated to this tree vertex

			//if the block does not have multiple important vertex and it does have important vertices --> it has exactly one important vertex!
			bool blockHasOnlyOneImportantVertex = !(bcTreeVertices[tv].flagHasMultipleImportantVertices()) && bcTreeVertices[tv].flagHasImportantVertex();
			
			if( blockHasOnlyOneImportantVertex &&   bcTreeVertices[tv].isLeaf()) { //is it a leaf with only one important vertex?
				//we should ignore if the only important vertex is an articulation STILL connecting to the important part of the tree...

				//we simply have to check if the important vertex is an articulation
				//and if the vertex of this articulation in the BT has been removed... 
				//(maybe it is an articulation, but an articulation that connects to the erased part of the tree
				// --- in this case the block is not ignored)
				bool shouldIgnore = false;
				
				
				int start = biconnectedCompsRaggedArrayStart[idComponent];
				int end = biconnectedCompsRaggedArrayStart[idComponent+1];
				for(int i=start;i<end;i++) { //for each vertex compV in this biconnected component
					int compV = biconnectedCompsRaggedArray[i]; 
					if(graph.isImportantVertex(compV)) {
						//is this important vertex an articulation and this articulation is still
						//in the tree? --> ignore component!
						if(isArticulation[compV] && !bcTreeVertices[blockId[compV]].flagDisabled())
							shouldIgnore = true;
						break;
					}
					
				} 
				if(shouldIgnore) continue;
			}
			//special case: 
			//we have a single vertex in this BC tree and component has only one vertex... (if this vertex is both starting pt and controller --> we have already selected it to be written (but others in the componnent shouldnt be written))
			if(bcTreeVertices[tv].numActiveNeighbors==0 && blockHasOnlyOneImportantVertex) continue;

			//all the vertices in this biconnected component are in a simple path between a controller and a starting point	
			int start = biconnectedCompsRaggedArrayStart[idComponent];
			int end = biconnectedCompsRaggedArrayStart[idComponent+1];
			outputVerticesIds.insert(outputVerticesIds.end(),biconnectedCompsRaggedArray.begin()+start,biconnectedCompsRaggedArray.begin()+end);
		}
	}

	//ensure we have only unique vertices ids in the output...
	{
			LOGGING(Timer t0("sort the ids (the internal ids used in our graph... i.e,. 0,1,2,...nv-1) output vertices and get the unique ones");)
			sort(outputVerticesIds.begin(),outputVerticesIds.end());
			outputVerticesIds.resize(unique(outputVerticesIds.begin(),outputVerticesIds.end())-outputVerticesIds.begin());
	}

}


