all:
	clear
	g++ -std=c++0x -fopenmp fchc.cc
	export OMP_NUM_THREADS=12
	./a.out
