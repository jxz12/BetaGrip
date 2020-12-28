randomkit.o: randomkit.c randomkit.h
	gcc randomkit.c -c

grip.o: grip.cpp grip.hpp
	g++ grip.cpp -std=c++11 -c

main.o: main.cpp grip.hpp
	g++ main.cpp -std=c++11 -c

main: randomkit.o grip.o main.o
	g++ randomkit.o grip.o main.o -o main

run: main ../christmas\ carol.txt
	./main ../christmas\ carol.txt

clean: 
	rm main *.o