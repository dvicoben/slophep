#basic makefile
CC=g++
#CC=g++ -std=c++11
CFLAGS=-c -Wall

#all: SL_Decay

BdToDstFFBGL: BdToDstFFBGL.o SL_Decay.o
	$(CC) BdToDstFFBGL.o SL_Decay.o -o BdToDstFFBGL

BdToDstFFBGL.o: BdToDstFFBGL.cpp
	$(CC) $(CFLAGS) BdToDstFFBGL.cpp

BsToDsstFFBGL: BsToDsstFFBGL.o SL_Decay.o
	$(CC) BsToDstFFBGL.o SL_Decay.o -o BsToDsstFFBGL

BsToDsstFFBGL.o: BsToDsstFFBGL.cpp
	$(CC) $(CFLAGS) BsToDsstFFBGL.cpp

SL_Decay: SL_DecayApp.o SL_Decay.o 
	$(CC) SL_DecayApp.o SL_Decay.o -o SL_Decay

SL_DecayApp.o: SL_DecayApp.cpp
	$(CC) $(CFLAGS) SL_DecayApp.cpp

SL_Decay.o: SL_Decay.h SL_Decay.cpp SL_Decay_Config.h SL_Decay_Config.cpp
	$(CC) $(CFLAGS) SL_Decay.cpp

clean:
	rm -rf *.o  *~ hellomake 