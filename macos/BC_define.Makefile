# ─────────────────────────────────────────────────────────────
#  Makefile – FlowState BC_define & helper tools
#  Apple-silicon / Homebrew CGNS + HDF5
# ─────────────────────────────────────────────────────────────

# Make "all" the default target so 'make' builds both tools.
.DEFAULT_GOAL := all

# ── compiler / flags ─────────────────────────────────────────
CXX      := clang++
CXXFLAGS := -O3 -march=native -ffast-math -std=c++20 \
            -Iinclude -I/opt/homebrew/include \
            -DVERIFY_TRANSFORM
LDFLAGS  := -L/opt/homebrew/lib -lcgns -lhdf5

# ── directory layout ────────────────────────────────────────
SRC_DIR  := src
BUILD_DIR:= build
BIN_DIR  := bin
SCRIPT_DIR := scripts
K_VALUES_H := include/k_values.hpp

K_SCRIPT := $(SCRIPT_DIR)/compute_k_values.sh

# Create build / bin automatically
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

# ── source files per executable ─────────────────────────────
# (add or remove *.cpp files here as the project grows)
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

# ── pattern rules ───────────────────────────────────────────
define make-exe
$1_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$($1_SRC))
$2 : $$($1_OBJS)
	$(CXX) $$^ $(LDFLAGS) -o $$@
endef

$(eval $(call make-exe,BC_DEFINE,$(BIN_DIR)/bc_define))
$(eval $(call make-exe,BC_DUMP,  $(BIN_DIR)/bc_dump))

# ── generic build rule ──────────────────────────────────────
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp $(K_VALUES_H)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ── phony targets ───────────────────────────────────────────
.PHONY: all clean regen_k_values
all: $(BIN_DIR)/bc_define $(BIN_DIR)/bc_dump

regen_k_values: $(K_SCRIPT)
	$(K_SCRIPT) > $(K_VALUES_H)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Cleaned."

# ─────────────────────────────────────────────────────────────
#  Usage:
#     make          – builds bc_define & bc_dump
#     make clean    – removes build/ and bin/ directories
# ─────────────────────────────────────────────────────────────
