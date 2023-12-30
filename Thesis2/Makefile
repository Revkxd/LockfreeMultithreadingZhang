CC = gcc
CFLAGS = -Wall

TARGET = seq
EXE_NAME = run_seq

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) $(TARGET).c -o $(EXE_NAME)

clean:
	rm -f $(EXE_NAME)