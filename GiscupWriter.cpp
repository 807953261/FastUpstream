
//Given the ids of the vertices that are upstream from a controller (i.e., in a path between a controller and a starting point)
//this function writes to the stream the edges and the vertices that are upstream from a controller
//if the two vertices of an edge e are in a path --> e will be in a path --> we write "e"
template <class GlobalIdType, class ParserType>
void writeGlobalIdsOfVerticesAndEdgesInPaths(ParserType &parser, const GiscupGraph &graph, vector<int> &outputVerticesIds, ostream &outputStream) {
	LOGGING(Timer t("Total writing, unique, saving disk, etc"));
	LOGGING(clog << "Number of output vertices ids: " << outputVerticesIds.size() << endl;)
	int nv = graph.numVertices();

	size_t numVertIdsWritten = 0;
	vector< int > outputEdgesIds; //ids of the edges that will be written...
	{
		LOGGING(Timer t("extract global ids, get edges global ids and write everything to disk...");)

		//can't use bool in this vector (race condition)
		vector<int> isInOutput(nv,0); //isInOutput[v] --> is graph vertex v in the output?
		int numOutputVerticesIds = outputVerticesIds.size();		
		
		#pragma omp parallel
		{
			ostringstream os; //each tread writes to a private ostringstream (so that the output is generated in parallel safely and w/o the need of synchronization)
			LOGGING(int myIdsWritten = 0;)

			#pragma omp for
			for(int i=0;i<numOutputVerticesIds;i++) {
				int idInOutput = outputVerticesIds[i];
				isInOutput[idInOutput] = true;

				//artificial vertices are not written to the output
				//if an artificial vertex is in the output --> at least one of its incident artificial edges (with same global id) will be (thus, we would have duplicate values)
				//the ids of this edge will be written (when we write the edges below...)
				if(!parser.isArtificialVertex(idInOutput)) {
					LOGGING(myIdsWritten++;)
					parser.getIdentifierVertex(idInOutput).print(os);
					os << "\n";					
				}
			}

			#pragma omp critical
			{
				outputStream << os.str();
				LOGGING(numVertIdsWritten += myIdsWritten;)
			}
		}
		

		vector<bool> isEdgeInOutput(graph.numEdges(),false);
		outputEdgesIds.reserve(numVertIdsWritten); //we don't know the amount of edges in the output yet... this is just an heuristic to (try to) save time
		for(int v:outputVerticesIds) {
			const int nv = graph.vertexDegree(v);
			for(int i=0;i<nv;i++) {
				int n = graph.getAdjVertex(v,i);
				//v,n is the i-th edge in v...
				if(v<=n) //we always have two copies of each edge: (u,v) and (v,u) ... let's write only one of them...
					if(isInOutput[n]) {
						//both sides of the i-th edge in v are in the output... write the edge!
						int edgeId = graph.getAdjEdge(v,i);
						if(!isEdgeInOutput[edgeId]) { //for uniqueness...
							isEdgeInOutput[edgeId] = true;
							outputEdgesIds.push_back(edgeId);
						}
					}
			}
		}
	}
	//now, write the edge output global ids in the file 
	{ 

		#pragma omp parallel
		{
			ostringstream os;

			#pragma omp for
			for(int i=0;i<outputEdgesIds.size();i++) {
				parser.getIdentifierEdge(outputEdgesIds[i]).print(os);
				os << "\n";	
			}

			#pragma omp critical
			{
				outputStream << os.str();
			}
		}
	}
	LOGGING(clog << "Number of globalIds written: " << outputEdgesIds.size() + numVertIdsWritten<< endl;)
}

