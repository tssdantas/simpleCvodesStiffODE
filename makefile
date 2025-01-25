CC = g++
CFLAGS = -I/usr/local/include
LDFLAGS = -L/usr/local/lib -lsundials_cvodes -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense

solverode: solverode.cpp
	$(CC) $(CFLAGS) SolverODE.cpp $(LDFLAGS) -o solverode

clean:
	rm -f solverode
