# ─────────────────────────────────────────────────────────────
#  Makefile – FlowState BC_define & helper tools (Linux)
# ─────────────────────────────────────────────────────────────

.DEFAULT_GOAL := all

CXX      ?= g++
CXXFLAGS := -O3 -march=native -ffast-math -std=c++20 \
            -Iinclude -I$(CGNS_INC) \
            -DVERIFY_TRANSFORM
LDFLAGS  := -L$(CGNS_LIB) -lcgns -lhdf5

SRC_DIR  := src
BUILD_DIR:= build
BIN_DIR  := bin
SCRIPT_DIR := scripts
K_VALUES_H := include/k_values.hpp

K_SCRIPT := $(SCRIPT_DIR)/compute_k_values.sh

$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

BC_DEFINE_SRC := $(SRC_DIR)/bc_define_main.cpp \
                 $(SRC_DIR)/bc_define_bcspec.cpp \
                 $(SRC_DIR)/bc_define_cgns.cpp \
                 $(SRC_DIR)/bc_define_plot3d.cpp \
                 $(SRC_DIR)/bc_define_helpers.cpp \
                 $(SRC_DIR)/mesh_io.cpp \
                 $(SRC_DIR)/cgns_cleanup.cpp \
                 $(SRC_DIR)/plot3d_io.cpp \
                 $(SRC_DIR)/connectivity_core.cpp \
                 $(SRC_DIR)/connectivity_face_build.cpp \
                 $(SRC_DIR)/connectivity_hash.cpp \
                 $(SRC_DIR)/connectivity_verify.cpp \
                 $(SRC_DIR)/connectivity_utils.cpp \
                 $(SRC_DIR)/connectivity_writer.cpp \
                 $(SRC_DIR)/logger.cpp

BC_DUMP_SRC   := $(SRC_DIR)/bc_dump.cpp \
                 $(SRC_DIR)/mesh_io.cpp \
                 $(SRC_DIR)/cgns_cleanup.cpp \
                 $(SRC_DIR)/logger.cpp

define make-exe
$1_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$($1_SRC))
$2 : $$($1_OBJS)
	$(CXX) $$^ $(LDFLAGS) -o $$@
endef

$(eval $(call make-exe,BC_DEFINE,$(BIN_DIR)/bc_define))
$(eval $(call make-exe,BC_DUMP,  $(BIN_DIR)/bc_dump))

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp $(K_VALUES_H)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: all clean regen_k_values
all: $(BIN_DIR)/bc_define $(BIN_DIR)/bc_dump

regen_k_values: $(K_SCRIPT)
	$(K_SCRIPT) > $(K_VALUES_H)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Cleaned."

# Usage:
#   make -C BC_define -f ../linux/BC_define.Makefile CGNS_INC=/path/include CGNS_LIB=/path/lib
