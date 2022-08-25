all: prova
prova: Main.o Planetesimal.o
	g++ Main.o Planetesimal.o -o prova
Main.o: Main.cpp Planetesimal.h
	g++ -c Main.cpp -o Main.o
Planetesimal.o : Planetesimal.cpp Planetesimal.h
	g++ -c Planetesimal.cpp -o Planetesimal.o

all: number
number: Number_Planetesimal.o Planetesimal.o
	g++ Number_Planetesimal.o Planetesimal.o -o number
Number_Planetesimal.o: Number_Planetesimal.cpp Planetesimal.h
	g++ -c Number_Planetesimal.cpp -o Number_Planetesimal.o


