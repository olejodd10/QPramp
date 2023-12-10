sources = [
    "../src/qpramp.c",
    "../src/vector.c",
    "../src/matrix.c",
    "../src/iterable_set.c",
    "../src/indexed_vectors.c",
];

% The "-g" flag disables the -O2 flag that will otherwise overwrite -Ofast
mex -setup C;
mex("-I../include", "-g", "CFLAGS=$CFLAGS -Ofast", sources{:}, "qpramp_solve_mpc.c", "-output", "qpramp_solve_mpc");
mex("-I../include", "-g", "CFLAGS=$CFLAGS -Ofast", sources{:}, "qpramp_solve.c", "-output", "qpramp_solve");
