TARGET = nBodySim
SRC_DIR = src
INC_DIR = include

CXX = h5c++
NVCC = nvcc
CUDA_ARCH_FLAGS = -arch=sm_75

CXXFLAGS = -O3 -std=c++17 -I$(INC_DIR)
NVCCFLAGS = -O3 $(CUDA_ARCH_FLAGS) --std=c++17 -I$(INC_DIR)
LDFLAGS = $(shell pkg-config --cflags --libs glew glfw3) -lGL -lcudart

ALL_CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
MAIN_CPP_SRCS = $(filter-out $(SRC_DIR)/datacreator.cpp, $(ALL_CPP_SRCS))
CU_SRCS = $(wildcard $(SRC_DIR)/*.cu)

OBJS_CPP = $(MAIN_CPP_SRCS:.cpp=.o)
OBJS_CU  = $(CU_SRCS:.cu=.o)
ALL_OBJS = $(OBJS_CPP) $(OBJS_CU)

all: $(TARGET)

$(TARGET): $(ALL_OBJS)
	$(CXX) $(ALL_OBJS) -o $@ $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)

run: all
	./$(TARGET)

.PHONY: all clean run