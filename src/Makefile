CC		= g++
CC_FLAGS	= -g3 -Wall -O3 -I/usr/local/include
LD_FLAGS	= -L/usr/local/lib -lgsl -lgslcblas
REC_OBJECTS	= FKHoptm.o emission.o common-functions.o PropMod.o minimization.o

rec: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o FKHoptm $(LD_FLAGS)

FKHoptm.o: FKHoptm.cpp
	$(CC) $(CC_FLAGS) -c FKHoptm.cpp

emission.o: emission.cpp
	$(CC) $(CC_FLAGS) -c emission.cpp

common-functions.o: common-functions.cpp
	$(CC) $(CC_FLAGS) -c common-functions.cpp

PropMod.o: PropMod.cpp
	$(CC) $(CC_FLAGS) -c PropMod.cpp

minimization.o: minimization.cpp
	$(CC) $(CC_FLAGS) -c minimization.cpp
