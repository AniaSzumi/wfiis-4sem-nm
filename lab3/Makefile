CC=g++
CFLAGS=  -I/usr/local/include/ -lgsl -lgslcblas -MMD -MP

OUT=main

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
DEP=$(SRC:.cpp=.d)

$(OUT): $(OBJ)
	$(CC)  -o $@  $^ $(CFLAGS)
	

	
.PHONY: clean	
clean:
	@rm -f $(OUT) $(OBJ) $(DEP)
	
run:
	./main
