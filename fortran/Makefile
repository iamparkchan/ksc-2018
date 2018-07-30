FC = mpif90
TARGET = gs
OBJECTS = main.o gram_schmidt.o param.o
FFLAGS = -O0 -fopenmp
LIBS = 

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

gram_schmidt.o: param.o
main.o: gram_schmidt.o param.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -rf *.mod $(TARGET) $(OBJECTS)
