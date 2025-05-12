CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -pthread
LDFLAGS = -lm
EXECUTABLE = newton
SRC = newton.c

all: $(EXECUTABLE)

$(EXECUTABLE): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLE) *.o

.PHONY: all clean
