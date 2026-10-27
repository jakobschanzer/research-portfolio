FC = gfortran
FFLAGS ?= -O2 -std=legacy -Wall -Wextra
PYTHON ?= python3
JULIA ?= julia
BUILD_DIR := build

.PHONY: check test figures clean privacy-check

check: $(BUILD_DIR)/janal $(BUILD_DIR)/hbb $(BUILD_DIR)/ring_angle privacy-check

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/janal: src/fortran/host_analysis/janal.f src/fortran/host_analysis/vol.f | $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $^

$(BUILD_DIR)/hbb: src/fortran/hydrogen_bonds/hb4.f | $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $<

$(BUILD_DIR)/ring_angle: src/fortran/ring_geometry/ring_angle.f | $(BUILD_DIR)
	$(FC) $(FFLAGS) -o $@ $<

privacy-check:
	$(PYTHON) scripts/check_publication_boundary.py

test:
	$(PYTHON) -m unittest discover -s tests -v

figures:
	$(JULIA) --project=. analysis/occupancy/analyze_drying_transition.jl
	$(JULIA) --project=. analysis/pmf_decomposition/make_figure.jl
	$(JULIA) --project=. analysis/spatial_maps/make_figure.jl
	$(JULIA) --project=. analysis/occupancy/make_figures.jl
	$(JULIA) --project=. analysis/temperature_response/make_figure.jl

clean:
	rm -rf $(BUILD_DIR)
