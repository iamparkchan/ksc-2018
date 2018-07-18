CC = mpicc
TARGET = gs
OBJECTS = main.o gram_schmidt.o
CFLAGS = -O0 -std=c99 -fopenmp 
LIBS = -lm

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean :
	rm -f $(TARGET) $(OBJECTS)
