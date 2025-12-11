TARGET        = nBodySim
REPLAY_TARGET = replay

SRC_DIR = src
INC_DIR = include

CXX  = h5c++
NVCC = nvcc

# ==== CUDA включаем по умолчанию (для коллег) ====
# На твоём компе без CUDA:  make USE_CUDA=0
USE_CUDA ?= 1

CXXFLAGS  = -O3 -std=c++17 -I$(INC_DIR)
LDFLAGS   = $(shell pkg-config --cflags --libs glew glfw3) -lGL

# ---------- CUDA-настройки (только если USE_CUDA=1) ----------
ifeq ($(USE_CUDA),1)
    CUDA_ARCH_FLAGS = -arch=sm_75
    NVCCFLAGS = -O3 $(CUDA_ARCH_FLAGS) --std=c++17 -I$(INC_DIR)

    # Макрос для кода, чтобы отличать сборку с CUDA
    CXXFLAGS += -DUSE_CUDA

    # Линкуем с CUDA runtime
    LDFLAGS += -lcudart

    CU_SRCS  = $(wildcard $(SRC_DIR)/*.cu)
    OBJS_CU  = $(CU_SRCS:.cu=.o)
else
    CU_SRCS  =
    OBJS_CU  =
endif

# ---------- Список всех .cpp ----------
ALL_CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# ---------- .cpp для основного симулятора (без datacreator.cpp и replay.cpp) ----------
SIM_CPP_SRCS = $(filter-out \
    $(SRC_DIR)/datacreator.cpp \
    $(SRC_DIR)/replay.cpp, \
    $(ALL_CPP_SRCS))

SIM_OBJS = $(SIM_CPP_SRCS:.cpp=.o)

# ---------- .cpp для реплея ----------
REPLAY_SRCS = \
    $(SRC_DIR)/replay.cpp \
    $(SRC_DIR)/object.cpp \
    $(SRC_DIR)/renderer.cpp \
    $(SRC_DIR)/control.cpp

REPLAY_OBJS = $(REPLAY_SRCS:.cpp=.o)

# ---------- Таргеты по умолчанию ----------
all: $(TARGET)

# ---------- Основной бинарник с симуляцией ----------
$(TARGET): $(SIM_OBJS) $(OBJS_CU)
	$(CXX) $(SIM_OBJS) $(OBJS_CU) -o $@ $(LDFLAGS)

# ---------- Бинарник реплея ----------
$(REPLAY_TARGET): $(REPLAY_OBJS)
	$(CXX) $(CXXFLAGS) $(REPLAY_OBJS) -o $@ $(LDFLAGS)

# ---------- Общие правила компиляции ----------
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifeq ($(USE_CUDA),1)
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET) $(REPLAY_TARGET)

run: all
	./$(TARGET)

.PHONY: all clean run $(REPLAY_TARGET)
