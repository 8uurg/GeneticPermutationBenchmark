CC = g++
gomea: GOMEA.c
ifeq ($(OS), Windows_NT)
	$(CC) -Wall -O3 -c -DBUILD_DLL GOMEA.c
	$(CC) -Wall -O3 -shared  GOMEA.o -o GOMEA.dll
else
	$(CC) -Wall -O3 -o GOMEA.so -fPIC -shared GOMEA.c
endif

clean:
ifeq ($(OS), Windows_NT)
	rm GOMEA.dll
	rm GOMEA.o
else
	rm GOMEA.so
endif
