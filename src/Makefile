TARGET = simplex
OS := $(shell uname)

ifeq ($(OS), Darwin)
CXX = clang++
CXXFLAGS = -Wno-format -O3
LDFLAGS = -v -framework Accelerate
endif

ifeq ($(OS), Linux)
CXX = g++
CXXFLAGS = -Wno-format -I./extern/include/
LDLAGS = -v -L./extern/lib/ -lblas -llapack -llapacke
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

