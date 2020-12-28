randomkit.o: randomkit.c
	gcc randomkit.c -c

grip.o: grip.cpp
	g++ grip.cpp -std=c++11 -c

main.o: main.cpp
	g++ main.cpp -std=c++11 -c

main: randomkit.o grip.o main.o
	g++ randomkit.o grip.o main.o -o main

run: main ../christmas\ carol.txt
	./main ../christmas\ carol.txt