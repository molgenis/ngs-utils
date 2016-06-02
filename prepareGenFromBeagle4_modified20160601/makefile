#compiler
CXX=g++

#compiler flags
COPTI= -O3
CDEBG= -g
CWARN= -Wall -Wextra -Wno-sign-compare
CBAMH= -D_LIBBAM
CSHV1=

#linker flags
LOPTI= -O3
LDEBG= -g
LSTDD= -lm -lpthread -lboost_iostreams -lboost_program_options
LSTDS= -Wl,-Bstatic -lboost_iostreams -lz -lbz2 -Wl,-Bdynamic -lm -lpthread
LCL3S= -Wl,-Bstatic -lz -lbz2 -Wl,-Bdynamic -lm -lpthread
LBAMD= -lbam -lz

#executable file
EFILE= bin/prepareGenFromBeagle4

#header files
HFILE= $(shell find src -name *.h)

#source files
CFILE= $(shell find src -name *.cpp)

#source path
VPATH= $(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#include path
ISTDP= -Isrc
IBAMP= -Ilib
#ICL3P= -I/apps/software/Boost/1.57.0-goolf-1.7.20-Python-2.7.9/include
#ICL3P= -I$EBROOTBOOST/include
ICL3P= -I/apps/software/Boost/1.60.0-foss-2015b-Python-2.7.9/include/

#library path
LSTDP= -Llib

#object files
OFILE= $(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
OBOST=

#default
all: cluster

#dynamic release
dynamic: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
dynamic: LFLAG=$(LOPTI) $(LSTDD)
dynamic: IFLAG=$(ISTDP)
dynamic: $(EFILE)

#debug release
debug: CFLAG=$(CDEBG) $(CWARN) $(CSHV1)
debug: LFLAG=$(LDEBG) $(LSTDD)
debug: IFLAG=$(ISTDP)
debug: $(EFILE)

#static release
static: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
static: LFLAG=$(LOPTI) $(LSTDS)
static: IFLAG=$(ISTDP)
static: $(EFILE)

#cluster release
cluster: CFLAG=$(COPTI) $(CWARN) $(CSHV1)
cluster: LFLAG=$(LOPTI) $(LCL3S)
cluster: IFLAG=$(ISTDP) $(ICL3P)
#cluster: OBOST=/apps/software/Boost/1.57.0-goolf-1.7.20-Python-2.7.9/lib/libboost_iostreams.a /apps/software/Boost/1.57.0-goolf-1.7.20-Python-2.7.9/lib/libboost_program_options.a
OBOST=-lz /apps/software/Boost/1.60.0-foss-2015b-Python-2.7.9/lib/libboost_iostreams.a
cluster: $(EFILE)

$(EFILE): $(OFILE)
	$(CXX) $^ $(OBOST) -o $@ $(LFLAG)

obj/%.o: %.cpp $(HFILE)
	$(CXX) -o $@ -c $< $(CFLAG) $(IFLAG)

clean: 
	rm -f obj/*.o $(EFILE)
	
test:
	cp $(EFILE) ~/$(EFILE).v10

oxford:
	cp $(EFILE) ~/bin/.
	
install:
	cp $(EFILE) ./tmp/
