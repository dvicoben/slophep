#basic makefile
CC=g++
#CC=g++ -std=c++11
CFLAGS=-c -Wall

all: SL_Decay

SL_Decay: SL_DecayApp.o SL_Decay.o 
	$(CC) SL_DecayApp.o SL_Decay.o -o SL_Decay

SL_DecayApp.o: SL_DecayApp.cpp
	$(CC) $(CFLAGS) SL_DecayApp.cpp

SL_Decay.o: SL_Decay.h SL_Decay.cpp SL_Decay_Config.h SL_Decay_Config.cpp
	$(CC) $(CFLAGS) SL_Decay.cpp

clean:
	rm -rf *.o  *~ hellomake 