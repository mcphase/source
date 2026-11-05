# Build the macOS-native parser and packer, collect files from the Inno Setup
# script, create the compressed payload, build a small uninstaller, embed both
# resources, and link the final Windows installer.

SHELL := /bin/sh
.DEFAULT_GOAL := all

# Repository layout:
#   ~/mcphas/version
#   ~/mcphas/innosetup_mac.iss
#   ~/mcphas/bin/iss/Makefile
#   ~/mcphas/Output/
ROOT         := $(HOME)/mcphas
ISS_DIR      := $(ROOT)/bin/iss
VERSION_FILE := $(ROOT)/version
ISS_SCRIPT   := $(ROOT)/innosetup_mac.iss
OUTPUT_DIR   := $(ROOT)/Output

# All intermediate files are kept here, so `make clean` can remove them
# without touching the final installer in Output/.
BUILD_DIR    := $(ISS_DIR)/build
MANIFEST     := $(BUILD_DIR)/files-manifest.tsv
PAYLOAD      := $(BUILD_DIR)/payload.mcp
UNINSTALLER  := $(BUILD_DIR)/uninstall_stub.exe
RESOURCE_RC  := $(BUILD_DIR)/installer_resources.rc
RESOURCE_OBJ := $(BUILD_DIR)/installer_resources.o

# Native compiler for the parser and packer; MinGW-w64 tools cross-compile
# Windows executables while running the build on macOS.
CXX       ?= c++
MINGW_CXX ?= x86_64-w64-mingw32-g++
WINDRES   ?= x86_64-w64-mingw32-windres
CC ?= clang
MINGW_CC ?= x86_64-w64-mingw32-gcc

WIN_MINIZ_OBJECTS := $(BUILD_DIR)/win_miniz.o \
                     $(BUILD_DIR)/win_miniz_tdef.o \
                     $(BUILD_DIR)/win_miniz_tinfl.o

$(BUILD_DIR)/win_miniz.o: $(ISS_DIR)/miniz.c $(ISS_DIR)/miniz.h | $(BUILD_DIR)
	$(MINGW_CC) -O2 -c "$<" -o "$@"

$(BUILD_DIR)/win_miniz_tdef.o: $(ISS_DIR)/miniz_tdef.c $(ISS_DIR)/miniz_tdef.h | $(BUILD_DIR)
	$(MINGW_CC) -O2 -c "$<" -o "$@"

$(BUILD_DIR)/win_miniz_tinfl.o: $(ISS_DIR)/miniz_tinfl.c $(ISS_DIR)/miniz_tinfl.h | $(BUILD_DIR)
	$(MINGW_CC) -O2 -c "$<" -o "$@"


# Miniz's C implementation is compiled separately as C, not as C++.
MINIZ_OBJECTS := $(BUILD_DIR)/miniz.o \
                 $(BUILD_DIR)/miniz_tdef.o \
                 $(BUILD_DIR)/miniz_tinfl.o

$(BUILD_DIR)/miniz.o: $(ISS_DIR)/miniz.c $(ISS_DIR)/miniz.h | $(BUILD_DIR)
	$(CC) -std=c99 -O2 -c "$<" -o "$@"

$(BUILD_DIR)/miniz_tdef.o: $(ISS_DIR)/miniz_tdef.c $(ISS_DIR)/miniz_tdef.h | $(BUILD_DIR)
	$(CC) -std=c99 -O2 -c "$<" -o "$@"

$(BUILD_DIR)/miniz_tinfl.o: $(ISS_DIR)/miniz_tinfl.c $(ISS_DIR)/miniz_tinfl.h | $(BUILD_DIR)
	$(CC) -std=c99 -O2 -c "$<" -o "$@"


# Link the Windows installer with the same miniz objects and embedded resources.
$(INSTALLER): $(ISS_DIR)/installer_stub.cpp $(ISS_DIR)/miniz.h \
              $(MINIZ_OBJECTS) $(RESOURCE_OBJ) | $(OUTPUT_DIR)
	$(MINGW_CXX) $(WIN_FLAGS) -o "$@" \
		"$(ISS_DIR)/installer_stub.cpp" $(MINIZ_OBJECTS) \
		"$(RESOURCE_OBJ)" -lshell32 -lole32 -luuid

# Extract the numeric version from a line such as:
#   #define MCPHASVERSION "mcphas 6.1"
VERSION := $(shell sed -nE \
	's/.*MCPHASVERSION[[:space:]]+"[^"]*([0-9]+\.[0-9]+)".*/\1/p' \
	"$(VERSION_FILE)" | head -n 1)

ifeq ($(strip $(VERSION)),)
$(error Could not extract a version such as 6.1 from $(VERSION_FILE))
endif

# Convert dots to underscores for the installer filename:
#   6.1 -> McPhaseSetup6_1.exe
VERSION_TAG := $(subst .,_,$(VERSION))
INSTALLER   := $(OUTPUT_DIR)/McPhaseSetup$(VERSION_TAG).exe

# Native C++17 build flags and static MinGW runtime flags.
NATIVE_FLAGS := -std=c++17 -O2 -Wall -Wextra
WIN_FLAGS    := -std=c++17 -O2 -static -static-libgcc -static-libstdc++ \
                -DMCPHASE_VERSION='"$(VERSION)"'

.PHONY: all clean info

# The default target builds the versioned Windows installer.
all: $(INSTALLER)

# Print the version and final output path without building anything.
info:
	@echo "Version:   $(VERSION)"
	@echo "Installer: $(INSTALLER)"

# Create output directories when needed.
$(BUILD_DIR) $(OUTPUT_DIR):
	mkdir -p "$@"

# Build the native Inno Setup script parser.
$(BUILD_DIR)/iss_files: $(ISS_DIR)/iss_files.cpp | $(BUILD_DIR)
	$(CXX) $(NATIVE_FLAGS) -o "$@" "$<"

# Build the native packer, linking the vendored miniz implementation.
$(BUILD_DIR)/package: $(ISS_DIR)/package.cpp $(ISS_DIR)/miniz.h $(MINIZ_OBJECTS)
	$(CXX) $(NATIVE_FLAGS) -o "$@" \
		"$(ISS_DIR)/package.cpp" $(MINIZ_OBJECTS)

# Parse the .iss script and collect its file entries into a TSV manifest.
$(MANIFEST): $(BUILD_DIR)/iss_files $(ISS_SCRIPT)
	"$(BUILD_DIR)/iss_files" "$(ISS_SCRIPT)" "$@"

# Compress the files listed in the manifest into the installer payload.
$(PAYLOAD): $(BUILD_DIR)/package $(MANIFEST)
	"$(BUILD_DIR)/package" "$(MANIFEST)" "$@"

# Build the standalone Windows uninstaller. It does not include the payload
# or miniz, so it remains much smaller than the installer.
$(UNINSTALLER): $(ISS_DIR)/uninstaller_stub.cpp | $(BUILD_DIR)
	$(MINGW_CXX) $(WIN_FLAGS) -o "$@" "$<" \
		-lshell32 -lole32 -luuid

# Create a Windows resource script embedding the payload as resource 101
# and the small uninstaller as resource 102.
$(RESOURCE_RC): $(PAYLOAD) $(UNINSTALLER) | $(BUILD_DIR)
	printf '101 RCDATA "%s"\n102 RCDATA "%s"\n' \
		"$(abspath $(PAYLOAD))" "$(abspath $(UNINSTALLER))" > "$@"

# Compile the resource script for linking into the Windows installer.
$(RESOURCE_OBJ): $(RESOURCE_RC)
	$(WINDRES) "$<" -O coff -o "$@"

# Link the installer stub, miniz, and embedded resources into the final
# versioned executable in ~/mcphas/Output/.
$(INSTALLER): $(ISS_DIR)/installer_stub.cpp $(ISS_DIR)/miniz.h \
              $(WIN_MINIZ_OBJECTS) $(RESOURCE_OBJ) | $(OUTPUT_DIR)
	$(MINGW_CXX) $(WIN_FLAGS) -o "$@" \
		"$(ISS_DIR)/installer_stub.cpp" $(WIN_MINIZ_OBJECTS) \
		"$(RESOURCE_OBJ)" -lshell32 -lole32 -luuid
	@echo "Successfully produced $(notdir $(INSTALLER)) in $(OUTPUT_DIR) . Type make clean to remove all temporary files."

# Remove only generated intermediates. The final installer in Output/ is
# deliberately preserved.
clean:
	rm -rf "$(BUILD_DIR)"
