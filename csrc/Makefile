# Specify compiler
CC=cl.exe

# Specify linker
LINK=link.exe
TARGET = prog


$(TARGET): OptBlock.o
	cl $^ -o $@ -lm

OptBlock.o: OptBlock.c
	cl -c $< -o $@

lib.a: lib1.o lib2.o
	ar rcs $@ $^

lib1.o: lib1.c lib1.h
	cl -c -o $@ $<

lib2.o: lib2.c lib2.h
	gcc -c -o $@ $<

clean:
	rm -f *.o *.a $(TARGET)