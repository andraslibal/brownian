alpha:
	# icc alpha.cpp -o alpha.out -fopenmp -O3 -xAVX -fp-model fast=2 -ljsmn -L/home/bajcsiadel/jsmn -std=c++11
	g++ alpha.cpp -o alpha.out -O3 -ffast-math -mavx -ljsmn -L/home/bajcsiadel/jsmn -std=c++11

clean:
	rm alpha.out result
