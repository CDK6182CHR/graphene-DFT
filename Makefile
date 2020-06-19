target=graphene.exe
source=$(wildcard *.cpp)
OBJS=$(source:%.cpp=%.o)

CXX=g++
CC=gcc
CXXFLAGS=-g -Wall -w -std=c++11 -lgsl

$(target):$(OBJS)
	$(CXX) -o $(target) $(OBJS)  $(CXXFLAGS)
%.o:%.cpp
	$(CXX) -c -o $@ $<  $(CXXFLAGS)

.PHONY:clean
clean:
	del *.o $(target)