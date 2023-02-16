TARGET = simplex
DEBUG = 0
OS := $(shell uname)

INCLUDES = -I../boost_1_81_0/
LIBS =

ifeq ($(OS), Darwin)
CXX = clang++
CXXFLAGS = -std=c++17 $(INCLUDES)
LDFLAGS = -v -framework Accelerate
endif

ifeq ($(OS), Linux)
CXX = g++
CXXFLAGS = -Wno-format -std=c++17
LDFLAGS = -v -lopenblas
endif

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -DDEBUG
else
CXXFLAGS += -O3
endif

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

-include $(DEPS)

%.d: %.cpp
	$(CXX) $(CXXFLAGS) $< -MM -MP -MT $(@:.d=.o) >$@

clean:
	$(RM) $(SRCS:.cpp=.o) $(SRCS:.cpp=.d)
	$(RM) $(TARGET)

print-%  : ; @echo $* = $($*)
