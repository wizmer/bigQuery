INCLUDE=-I${HOME}/utils/include
AMSSTATICLIBS := -L$(AMSWD)/lib/linuxx8664gcc5.34 -lntuple_slc6_PG
LIBS=`root-config --cflags --libs ` -L${HOME}/utils/lib -lGeneralUtils
OPTIMIZATION=-g
mcmc:	mcmc.cpp pd_model.o
	g++ ${OPTIMIZATION} ${INCLUDE} mcmc.cpp -o mcmc -std=c++11 ${LIBS} pd_model.o

read:
	g++ ${OPTIMIZATION} ${INCLUDE} reader.cpp -o reader  ${LIBS} -lRootUtils $(AMSSTATICLIBS) 

pd_model.o: pd_model.cpp  pd_model.hpp
	g++ ${OPTIMIZATION} -c pd_model.cpp -std=c++11 


clean:
	rm -f mcmc reader pd_model.o

