randomkit.o: randomkit.c randomkit.h
	gcc randomkit.c -c

peach.o: peach.cpp peach.hpp
	g++ peach.cpp -std=c++11 -c

main.o: main.cpp peach.hpp
	g++ main.cpp -std=c++11 -c

main: randomkit.o peach.o main.o
	g++ randomkit.o peach.o main.o -o main

run: main ../christmas\ carol.txt
	./main ../christmas\ carol.txt

clean: 
	rm main *.o