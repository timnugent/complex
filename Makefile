CPP = g++
CPPFLAGS = -Wall -O3
INC = `pkg-config libgvc --cflags`
LIBS = -lboost_regex -lboost_graph -lboost_system -lboost_filesystem `pkg-config libcgraph libgvc --libs`

all: bin/gengraph

src/smithwaterman.o: src/smithwaterman.cpp 
	$(CPP) $(CPPFLAGS) -c src/smithwaterman.cpp -o src/smithwaterman.o
src/gengraph.o: src/gengraph.cpp
	$(CPP) $(INC) $(CPPFLAGS) -c src/gengraph.cpp -o src/gengraph.o
bin/gengraph: src/gengraph.o src/smithwaterman.o 
	$(CPP) $(CPPFLAGS) src/smithwaterman.o src/gengraph.o -o bin/gengraph ${LIBS} 
clean:
	rm bin/gengraph src/*.o examples/*.svg
test: examples/1kqf.pdb examples/1l0l.pdb examples/3pcq.pdb
	bin/gengraph -s 1 -l 2 -a B:1,C:4,D:1,F:4,H:1,I:4 -b C:8,D:1,E:1,G:1,J:1,K:1,N:8,O:1,P:1,R:1,U:1,V:1 examples/1kqf.pdb examples/1l0l.pdb
	bin/gengraph -s 1 -l 2 -t 1 -a A:11,B:11,F:1,G:1,H:1,I:2,J:3,K:1,L:1,M:11,N:11,R:1,S:1,T:1,U:2,V:3,W:1,Y:1,X:11,Z:11,4:1,5:1,6:1,7:2,8:3,9:1,a:1 examples/3pcq.pdb


