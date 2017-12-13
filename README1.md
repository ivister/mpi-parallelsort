# mpi-parallelsort
mpich
compile: mpic++ -o ParallelSort parallelsort.cpp -std=c++11 -Wall
run params: -f <filename> | -n <count of notes>
must be <count of notes> % (4 * np) == 0
