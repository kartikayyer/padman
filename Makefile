objects = recon.o excom.o max.o
CFLAGS = -c -Wall -O3

recon : $(objects)
	g++ -o recon $(objects) -Wall -O3
recon.o : cpprecon.cpp structs.hpp
	g++ -c cpprecon.cpp -o recon.o $(CFLAGS)
excom.o : excom.cpp structs.hpp
	g++ -c excom.cpp -o excom.o -lm $(CFLAGS)
max.o : max.cpp structs.hpp
	g++ -c max.cpp -o max.o -lm $(CFLAGS)
