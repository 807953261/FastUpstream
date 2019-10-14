flag= -O3 -std=c++11  -fopenmp -D_GLIBCXX_PARALLEL -DNDEBUG 


giscup.out: main.o GiscupSolver.o utils.o
	g++ $(flag)  main.o GiscupSolver.o utils.o -o  giscup.out

GiscupSolver.o: GiscupSolver.cpp GiscupSolver.h GiscupGraph.h logging.h GiscupGraph.h
	g++ $(flag) -c GiscupSolver.cpp   

utils.o: utils.h utils.cpp logging.h
	g++ $(flag) -c utils.cpp   

main.o: main.cpp GiscupWriter.cpp  GiscupGraph.h GiscupParser.h GiscupSolver.h  ImportantIdentifier.h RapidJsonParser.h RapidJsonParser.cpp logging.h
	g++ $(flag) -c main.cpp  -I libs/ 

