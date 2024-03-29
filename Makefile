CXX		  := g++
CXX_FLAGS := -Wall -O0 -g
RUN_ARGS  := --help

BIN		:= bin
BUILD   := build
SRC		:= src
INCLUDE	:= include
LIB		:= lib

LIBRARIES	:= -lboost_program_options
EXECUTABLE	:= merlin

SRCS := $(shell find $(SRC) -name *.cpp)
OBJS := $(SRCS:%=$(BUILD)/%.o)
DEPS := $(OBJS:.o=.d)

INC_FLAGS := $(addprefix -I,$(INCLUDE))
LDFLAGS := $(addprefix -L,$(LIB))

CPP_FLAGS := $(INC_FLAGS) -MMD -MP


$(BIN)/$(EXECUTABLE): $(OBJS)
	$(MKDIR_P) $(BIN)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@ $(LIBRARIES)

# c++ source
$(BUILD)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPP_FLAGS) $(CXX_FLAGS) -c $< -o $@


.PHONY: clean

all: $(BIN)/$(EXECUTABLE)

run: clean all
	./$(BIN)/$(EXECUTABLE) $(RUN_ARGS)

clean:
	$(RM) -r $(BUILD)/*
	$(RM) -r $(BIN)/*

-include $(DEPS)

MKDIR_P ?= mkdir -p

