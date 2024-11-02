CC = gcc
OPTIMIZE = -O3
DEBUG = -g
AR = ar rcs

CFLAGS = $(OPTIMIZE) $(DEBUG)

all: libfft.a

libfft.a: fft.o
	$(AR) $@ $^

fft.o: src/fft.c src/fft.h
	$(CC) -c $< -o $@ $(CFLAGS)


testfft: test/test.c libfft.a
	$(CC) $< -o $@ $(CFLAGS) -L . -lfft
	


.PHONY: clean

clean:
	$(RM) fft.o libfft.a *.exe
