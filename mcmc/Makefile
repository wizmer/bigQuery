INCLUDE=-I./utils
AMSSTATICLIBS := -L$(AMSWD)/lib/linuxx8664gcc5.34 -lntuple_slc6_PG
LIBS=`root-config --cflags --libs ` 
OPTIMIZATION=-O3

mcmc:	mcmc.cpp pd_model.o generalUtils.o
	g++ ${OPTIMIZATION} ${INCLUDE} $^ -o mcmc -std=c++11 ${LIBS}

read:
	g++ ${OPTIMIZATION} ${INCLUDE} reader.cpp -o reader  ${LIBS} -lRootUtils $(AMSSTATICLIBS) 

generalUtils.o: utils/generalUtils.cpp utils/generalUtils.hpp 
	g++ ${OPTIMIZATION} ${INCLUDE} -c $< -std=c++11 

pd_model.o: pd_model.cpp  pd_model.hpp
	g++ ${OPTIMIZATION} ${INCLUDE} -c $< -std=c++11 

test_matrix: tests/test_matrix.cpp utils/Matrix.hpp 
	g++ ${INCLUDE} -o $@ $< -std=c++11 

test_model: tests/test_model.cpp pd_model.o
	g++ ${INCLUDE} -o $@ $^ -std=c++11 

clean:
	rm -f mcmc reader pd_model.o

