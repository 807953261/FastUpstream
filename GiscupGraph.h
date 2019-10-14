#ifndef GISCUP_GRAPH_H
#define GISCUP_GRAPH_H

#include <cassert>
#include <vector>
#include <string>
#include <map>


struct Vertex {
	bool isController;
	bool isStartingPoint;
	int degree;
	int posAdjOfVStart; //this is used in the ragged array representation. The neighbors of vertex v will be stored after (including) position v.posAdjOfVStart
};

//This class basically represents the GISCUP graph
//It stores the adjacency information and some extra information about the vertices/edges (is starting point? is controller?)
//It is implemented as a ragged array (for performance)
class GiscupGraph {
public:
	GiscupGraph() {}

	//The degrees of the vertices must be known in advance (we use a ragged array...)
	//We assume the first "degrees.size()" vertices will have degrees as stored in degrees
	//the other ones will have degrees defaultDegreesOtherVertices
	GiscupGraph(int nv,const vector<size_t> &degrees, int defaultDegreesOtherVertices) {
		vertexInfo.resize(nv);

		allocateRaggedArray(degrees,defaultDegreesOtherVertices);
	}




	void addEdge(int u,int v,int edgeId) {
		assert(u<numVertices());
		assert(v<numVertices());
		int posInsert;

		//this allows adding edges in parallel safely
		#pragma omp atomic capture
		posInsert = vertexInfo[u].degree++;

		posInsert += vertexInfo[u].posAdjOfVStart;

		edges[posInsert][0] = v; 
		edges[posInsert][1] = edgeId; 
	}

	//an important vertex is a vertex that is a controller or a starting point
	bool isImportantVertex(int v) const { return vertexInfo[v].isController || isStartingPointVertexOrArtVertex(v);}
	bool isControllerVertex(int v) const { return vertexInfo[v].isController;}
	//is v the id of a starting point? (also returns true if v is the id of an artificial vertex representing a starting point that is an edge in the original input)
	bool isStartingPointVertexOrArtVertex(int v) const { return v<numVertices() && vertexInfo[v].isStartingPoint;}
	
	int numVertices() const { return vertexInfo.size(); }
	int vertexDegree(int v) const { return vertexInfo[v].degree; } //number of vertices adjacent to vertex v...
	int getAdjVertex(int v, int i) const { 
		assert(v<edges.size()); assert(i<vertexDegree(v)); 
		assert(vertexInfo[v].posAdjOfVStart+i < edges.size());
		return edges[vertexInfo[v].posAdjOfVStart+i][0]; 
	} //returns the id of the i-th neighbor of v...
	int getAdjEdge(int v, int i) const { return edges[vertexInfo[v].posAdjOfVStart+i][1];  } //get the id of the i-th edge of v


	void setStartingPointVertex(int vertexId, bool value) {
		assert(vertexId<numVertices());
		vertexInfo[vertexId].isStartingPoint = value;
	}
	void setControllerVertex(int vertexId, bool value) {
		assert(vertexId<numVertices());
		vertexInfo[vertexId].isController = value;
	}
	int numEdges() const { return edges.size(); }

private:	
	//each vertex is numbered from 0 to nv-1 (nv is the number of vertices)

	//edges stores a ragged array of edges
	std::vector<array<int,2> > edges;
	std::vector<Vertex> vertexInfo;


	//Allocates the memory for storing the edges
	//We assume the first "degrees.size()" vertices will have degrees as stored in degrees
	//the other ones will have degrees defaultDegreesOtherVertices
	void allocateRaggedArray(const vector<size_t> &degrees, int defaultDegreesOtherVertices) {		
		size_t numEdgesSeemSoFar = 0;
		int numDegrees= degrees.size();
		int nv = numVertices();

		for(int i=0;i<nv;i++) {
			vertexInfo[i].isController = vertexInfo[i].isStartingPoint = false;
			vertexInfo[i].posAdjOfVStart = numEdgesSeemSoFar;
			numEdgesSeemSoFar += (i<numDegrees)?degrees[i]:defaultDegreesOtherVertices;
		}
		edges.resize(0);
		edges.resize(numEdgesSeemSoFar);
	}
};

#endif