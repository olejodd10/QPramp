sources = [
    "../src/qp_ramp.c",
    "../src/vector.c",
    "../src/matrix.c",
    "../src/iterable_set.c",
    "../src/indexed_vectors.c",
];

% The "-g" flag disables the -O2 flag that will otherwise overwrite -Ofast
mex -setup C;
mex("-I../include", "-g", "CFLAGS=$CFLAGS -Ofast", sources{:}, "qp_ramp_solve_mpc.c", "-output", "qp_ramp_solve_mpc");
mex("-I../include", "-g", "CFLAGS=$CFLAGS -Ofast", sources{:}, "qp_ramp_solve.c", "-output", "qp_ramp_solve");