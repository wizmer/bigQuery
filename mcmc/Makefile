INCLUDE=-I./utils -I$(amsUtils)/include -I.
AMSSTATICLIBS := -L$(AMSWD)/lib/linuxx8664gcc5.34 -lntuple_slc6_PG
LIBS=`root-config --cflags --libs | sed 's/-lPhysics//g' | sed 's/-lGraf3d//g' | sed 's/-lPostscript//g' `  -L$(amsUtils)/lib
#LIBS=`root-config --cflags --libs`
OPTIMIZATION=-g
CXXFLAGS=-Wall -Weffc++ -pedantic -Wno-sign-compare \
    -pedantic-errors -Wextra -Wcast-align \
    -Wcast-qual  -Wchar-subscripts  -Wcomment \
    -Wdisabled-optimization \
    -Werror -Wfloat-equal  -Wformat  -Wformat=2 \
    -Wformat-nonliteral -Wformat-security  \
    -Wformat-y2k \
    -Wimport  -Winit-self  -Winline \
    -Winvalid-pch   \
    # -Wunsafe-loop-optimizations
-Wlong-long -Wmissing-braces \
    -Wmissing-field-initializers -Wmissing-format-attribute   \
    -Wmissing-include-dirs -Wmissing-noreturn \
    -Wparentheses  -Wpointer-arith \
    -Wredundant-decls -Wreturn-type \
    -Wsequence-point  -Wshadow -Wstack-protector \
    -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default \
    -Wswitch-enum -Wtrigraphs  -Wuninitialized \
    -Wunknown-pragmas  -Wunreachable-code -Wunused \
    -Wunused-function  -Wunused-label  -Wunused-parameter \
    -Wunused-value  -Wunused-variable  -Wvariadic-macros \
    -Wvolatile-register-var  -Wwrite-strings

all: mcmc test_model test_matrix

mcmc:	mcmc.cpp pd_model.o generalUtils.o
	g++ ${OPTIMIZATION} $(CXXFLAGS) ${INCLUDE} $(LIBS) $^ -o mcmc -std=c++11

reader:	reader.cpp
	g++ ${OPTIMIZATION} -std=c++11  -L${amsUtils}/lib ${INCLUDE} reader.cpp -o reader  ${LIBS} -lRootUtils  -lStack -lGeneralUtils

generalUtils.o: utils/generalUtils.cpp utils/generalUtils.hpp 
	g++ ${OPTIMIZATION} ${INCLUDE} -c $< -std=c++11 

pd_model.o: pd_model.cpp  pd_model.hpp
	g++ ${OPTIMIZATION} $(CXXFLAGS) ${INCLUDE} -c $< -std=c++11 

test_matrix: tests/test_matrix.cpp utils/Matrix.hpp  generalUtils.o
	g++ ${INCLUDE} $(CXXFLAGS) -o $@ $< -std=c++11 

test_model: tests/test_model.cpp pd_model.o generalUtils.o
	g++ ${INCLUDE} $(CXXFLAGS) -o $@ $^ -std=c++11 

clean:
	rm -f mcmc reader pd_model.o test_model

