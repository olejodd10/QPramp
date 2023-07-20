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

int assert_size_csv(FILE* f, size_t m, size_t n) {
    if (m != parse_matrix_csv_height(f) || n != parse_matrix_csv_width(f)) {
        return -1;
    }
    return 0;
}

static int parse_row_csv(FILE* f, size_t n, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        int ret = fscanf(f, "%lf,", &res[i]);
        if (ret != 1) {
            return -1;
        }
    }
    return 0;
}

int parse_vector_csv(FILE* f, size_t n, double res[n]) {
    int ret = assert_size_csv(f, n, 1); // CSVs with one row could just as well be accepted
    if (ret < 0) {
        return ret;
    }
    return parse_row_csv(f, n, res);
}

int parse_matrix_csv(FILE* f, size_t m, size_t n, double res[m][n]) {
    int ret = assert_size_csv(f, m, n);
    if (ret < 0) {
        return ret;
    }
    for (size_t i = 0; i < m; ++i) {
        ret = parse_row_csv(f, n, res[i]);
        if (ret < 0) {
            return ret;
        }
    }
    return 0;
}

