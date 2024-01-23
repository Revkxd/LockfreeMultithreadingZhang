CC = gcc
CFLAGS = 

SEQ = seq
HTSEQ = htseq

all: seq htseq

seq: $(SEQ).c
	$(CC) $(CFLAGS) $(SEQ).c -o $(SEQ)

htseq: $(HTSEQ).c
	$(CC) $(CFLAGS) $(HTSEQ).c -o $(HTSEQ)

clean:
	rm -f *.exe $(SEQ) $(HTSEQ)