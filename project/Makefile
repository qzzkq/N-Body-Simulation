TARGET        := nBodySim
REPLAY_TARGET := replay
CREATOR_TARGET := datacreator

SRC_DIR   := src
INC_DIR   := include
ALGO_DIR  := algo
BUILD_DIR := build

CXX  := h5c++
NVCC := nvcc

USE_CUDA ?= 1
CUDA_ARCH ?= -arch=sm_75

OPT_FLAGS := -O3 -march=native -ffast-math -flto

PKG_LIBS  := $(shell pkg-config --cflags --libs glew glfw3)
DEPFLAGS  := -MMD -MP

CXXFLAGS  := $(OPT_FLAGS) -std=c++17 -I$(INC_DIR) -I$(ALGO_DIR) $(DEPFLAGS)
LIBS      := $(PKG_LIBS) -lGL

ifeq ($(USE_CUDA),1)
    NVCCFLAGS := -O3 $(CUDA_ARCH) --std=c++17 -I$(INC_DIR) -I$(ALGO_DIR) -MD
    CXXFLAGS += -DUSE_CUDA
    LIBS     += -lcudart

    CU_SRCS  := $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(ALGO_DIR)/*.cu)
    OBJS_CU  := $(CU_SRCS:%.cu=$(BUILD_DIR)/%.o)
else
    CU_SRCS  :=
    OBJS_CU  :=
endif

ALL_CPP_SRCS := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(ALGO_DIR)/*.cpp)

SIM_CPP_SRCS := $(filter-out $(SRC_DIR)/datacreator.cpp $(SRC_DIR)/replay.cpp, $(ALL_CPP_SRCS))
SIM_OBJS     := $(SIM_CPP_SRCS:%.cpp=$(BUILD_DIR)/%.o)

REPLAY_SRCS  := $(SRC_DIR)/replay.cpp $(SRC_DIR)/object.cpp $(SRC_DIR)/renderer.cpp $(SRC_DIR)/control.cpp $(SRC_DIR)/physics.cpp
REPLAY_OBJS  := $(REPLAY_SRCS:%.cpp=$(BUILD_DIR)/%.o)

CREATOR_SRCS := $(SRC_DIR)/datacreator.cpp $(SRC_DIR)/object.cpp 
CREATOR_OBJS := $(CREATOR_SRCS:%.cpp=$(BUILD_DIR)/%.o)

DEPS := $(SIM_OBJS:.o=.d) $(REPLAY_OBJS:.o=.d) $(CREATOR_OBJS:.o=.d) $(OBJS_CU:.o=.d)

.PHONY: all clean run

all: $(TARGET) $(REPLAY_TARGET) $(CREATOR_TARGET)

$(TARGET): $(SIM_OBJS) $(OBJS_CU)
	$(CXX) $(CXXFLAGS) $(SIM_OBJS) $(OBJS_CU) -o $@ $(LIBS)

$(REPLAY_TARGET): $(REPLAY_OBJS)
	$(CXX) $(CXXFLAGS) $(REPLAY_OBJS) -o $@ $(LIBS)

$(CREATOR_TARGET): $(CREATOR_OBJS)
	$(CXX) $(CXXFLAGS) $(CREATOR_OBJS) -o $@ $(LIBS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifeq ($(USE_CUDA),1)
$(BUILD_DIR)/%.o: %.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif

-include $(DEPS)

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET) $(REPLAY_TARGET) $(CREATOR_TARGET)

run: all
	./$(TARGET)