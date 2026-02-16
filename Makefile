.PHONY: configure build run run-standing run-radial run-radial-absorbing mesh clean-solution

BUILD_DIR := cmake-build-release
BINARY := $(BUILD_DIR)/wave-equation
CONFIG ?= configs/standing.cfg

configure:
	cmake -S . -B $(BUILD_DIR)

build:
	@if [ ! -f $(BUILD_DIR)/CMakeCache.txt ]; then $(MAKE) configure; fi
	cmake --build $(BUILD_DIR) -j

run: build
	$(BINARY) --config $(CONFIG)

run-standing: build
	$(BINARY) --config configs/standing.cfg

run-radial: build
	$(BINARY) --config configs/radial.cfg

run-radial-absorbing: build
	$(BINARY) --config configs/radial_absorbing.cfg

mesh:
	./scripts/generate_mesh.sh

clean-solution:
	rm -f solution/*.vtu
