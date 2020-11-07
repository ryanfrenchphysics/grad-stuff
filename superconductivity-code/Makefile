C++FLAG = -lstdc++
CFLAGS = -lm -std=c++14 -I /usr/include/python3.7m/ -L/usr/local/lib -lpython3.7m -Wall -fconcepts -g
GCC = gcc
SRC = test.cpp Grid.cpp
OBJ = $(SRC:.cpp = .o)

Test: $(OBJ)
	$(GCC) $(C++FLAG) $(CFLAGS) -o Test $(OBJ)
clean:
	rm -f core *.o
