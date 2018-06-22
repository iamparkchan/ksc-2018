FC = mpif90
TARGET = gs
OBJECTS = main.o gram_schmidt.o
FFLAGS = -O0 -fopenmp
LIBS = 

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

main.o: gram_schmidt.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -rf *.mod $(TARGET) $(OBJECTS)
