CC = g++ 
all : boruvka prim kruskal

boruvka: boruvka.cpp
	$(CC) -o boruvka boruvka.cpp
	
prim: prim.cpp
	$(CC) -o prim prim.cpp
	
kruskal: kruskal.cpp
	$(CC) -o kruskal kruskal.cpp

clean:
	rm -f prim *~
	rm -f boruvka *~
	rm -f kruskal *~
