SIMULATION_TIMESTEPS := 60000
EPS := 2.2204e-15

EXAMPLE := 3
ifeq ($(EXAMPLE), 1)
    N := 2
    M := 1
    HORIZON := 10
    C := 66
else ifeq ($(EXAMPLE), 2)
    N := 4
    M := 2
    HORIZON := 30
    C := 316
else ifeq ($(EXAMPLE), 3)
    N := 25
    M := 4
    HORIZON := 30
    C := 1560
endif

INPUT_DIR_RAW := ../examples/example$(EXAMPLE)
INPUT_DIR := \"$(INPUT_DIR_RAW)\"

INPUT_VECTORS := x0 w
INPUT_MATRICES := a b invh g s f

WRITE_FORMAT := \"%.4le\"

DEFINES := SIMULATION_TIMESTEPS EPS N M HORIZON C INPUT_DIR WRITE_FORMAT
