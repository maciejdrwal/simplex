TARGET = simplex
DEBUG = 0
OS := $(shell uname)
RM = rm -f

INCLUDES = -Iinclude/
LIBS =

ifeq ($(OS), Darwin)
CXX = clang++
CXXFLAGS = -std=c++17
LDFLAGS = -v -framework Accelerate
endif

ifeq ($(OS), Linux)
CXX = g++
EXTRA_CXXFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0
CXXFLAGS = -std=c++17 $(EXTRA_CXXFLAGS) 
LDFLAGS = -ldl -lstdc++fs
endif

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -DDEBUG
else
CXXFLAGS += -O3 -flto
endif

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS) $(LDFLAGS) -Wl

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	$(RM) $(OBJS)
	$(RM) $(TARGET)

print-%  : ; @echo $* = $($*)
