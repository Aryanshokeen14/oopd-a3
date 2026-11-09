CXX = clang++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

TARGET = pandemic_sim
SRC = pandemic_sim.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET) *.csv

.PHONY: all clean

