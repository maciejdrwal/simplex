TARGET = simplex
DEBUG = 1
OS := $(shell uname)

ifeq ($(OS), Darwin)
CXX = clang++
CXXFLAGS = -Wno-format -std=c++14
LDFLAGS = -v -framework Accelerate
endif

ifeq ($(OS), Linux)
CXX = g++
CXXFLAGS = -Wno-format -std=c++14 -I../extern/include/
LDFLAGS = -v -L../extern/lib/ -lopenblas
endif

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -DDEBUG
else
CXXFLAGS += -O2
endif


src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
dep = $(obj:.o=.d)

$(TARGET) : $(obj)
	$(CXX) -o $@ $^ $(LDFLAGS)

-include $(dep)

$.d : $.cpp
	@$(CXX) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@

.PHONY: clean
clean:
	rm -rf $(obj) $(TARGET)

.PHONY: cleandep
cleandep:
	rm -rf $(dep)

print-%  : ; @echo $* = $($*)
