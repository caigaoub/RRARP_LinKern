CPP       	= g++
CPPARGS   	= -O3 -m64 -std=c++1y -Wall -Wextra -pedantic
SRCPATH		= ./src/
BINPATH		= ./bin/
DATPATH		= ./dat/
INCBOOST  	= -I /usr/include/
CPPLIBBOOST 	= -L /usr/lib/  -lboost_system -lboost_filesystem

all:
	$(CPP) $(CPPARGS) $(SRCPATH)DataHandler.cpp $(SRCPATH)PartitionScheme.cpp \
	-o $(BINPATH)main $(SRCPATH)main.cpp $(INCBOOST) $(CPPLIBBOOST)

clean:
	rm -rf $(BINPATH)*.o $(BINPATH)*.dSYM $(BINPATH)main

run:
	$(BINPATH)main
