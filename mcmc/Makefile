all: read mcmc

INCLUDE=-I${HOME}/utils/include
AMSSTATICLIBS := -L$(AMSWD)/lib/linuxx8664gcc5.34 -lntuple_slc6_PG
LIBS=`root-config --cflags --libs ` -L${HOME}/utils/lib -lRootUtils 
mcmc:
	g++ -g ${INCLUDE} mcmc.cpp -o main -std=c++11

read:
	g++ -g ${INCLUDE} reader.cpp -o reader  ${LIBS} $(AMSSTATICLIBS)

pd_model.o: pd_model.cpp  pd_model.hpp
	g++ -c pd_model.cpp -std=c++11 


clean:
	rm -f main reader

