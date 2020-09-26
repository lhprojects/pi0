FLAGS=$(shell root-config --libs) $(shell root-config --cflags)
pi0: pi0.cxx
	g++ pi0.cxx $(FLAGS)  blossom5-v2.05/blossom5.so  -o pi0
