#include "parse.h"

size_t parse_matrix_csv_height(FILE* f) {
    size_t count = 0;
    fpos_t pos;
    fgetpos(f, &pos);
    rewind(f);
    while (fscanf(f, "%*s") != EOF) {
        ++count;
    }
    fsetpos(f, &pos);
    return count;
}

size_t parse_matrix_csv_width(FILE* f) {
    size_t count = 0;
    fpos_t pos;
    fgetpos(f, &pos);
    rewind(f);
    while (fscanf(f, "%*f,") != EOF) {
        ++count;
        int c = fgetc(f);
        if (c == '\n' || c =='\r') {
            break;
        }
        ungetc(c, f);
    }
    fsetpos(f, &pos);
    return count;
}

void parse_vector_csv(FILE* f, size_t n, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        fscanf(f, "%lf,", &res[i]);
    }
}

void parse_matrix_csv(FILE* f, size_t m, size_t n, double res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        parse_vector_csv(f, n, res[i]);
    }
}

