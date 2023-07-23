#include "csv.h"

static size_t parse_matrix_csv_height(FILE* f) {
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

static size_t parse_matrix_csv_width(FILE* f) {
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

static int assert_size_csv(FILE* f, size_t m, size_t n) {
    if (m != parse_matrix_csv_height(f) || n != parse_matrix_csv_width(f)) {
        return -1;
    }
    return 0;
}

// Accepts both row and column vectors
static int parse_row_csv(FILE* f, size_t n, double res[n]) {
    for (size_t i = 0; i < n; ++i) {
        int ret = fscanf(f, "%lf,", &res[i]);
        if (ret != 1) {
            return -1;
        }
    }
    return 0;
}

int parse_vector_csv(const char* path, size_t n, double res[n]) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
    int ret = assert_size_csv(f, n, 1); // CSVs with one row could just as well be accepted
    if (ret < 0) {
        return ret;
    }
    ret = parse_row_csv(f, n, res);
    fclose(f);
    return ret;
}

int parse_matrix_csv(const char* path, size_t m, size_t n, double res[m][n]) {
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }
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
    fclose(f);
    return ret;
}

static ssize_t write_row_csv(FILE* f, size_t n, const double vec[n]) {
    size_t sum = 0;
    int ret = 0;
    for (size_t i = 0; i < n-1; ++i) {
        ret = fprintf(f, WRITE_FORMAT ",", vec[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = fprintf(f, WRITE_FORMAT "\n", vec[n-1]);
    if (ret < 0) {
        return ret;
    } else {
        sum += ret;
    }
    return sum;
}

// Writes as a column vector
ssize_t save_vector_csv(const char* path, size_t n, const double vec[n]) {
    FILE* f = fopen(path, "w+");
    if (f == NULL) {
        return -1;
    }
    size_t sum = 0;
    int ret;
    for (size_t i = 0; i < n; ++i) {
        ret = fprintf(f, WRITE_FORMAT "\n", vec[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = assert_size_csv(f, n, 1);
    if (ret < 0) {
        return ret;
    }
    fclose(f);
    return sum;
}

ssize_t save_matrix_csv(const char* path, size_t m, size_t n, const double mat[m][n]) {
    FILE* f = fopen(path, "w+");
    if (f == NULL) {
        return -1;
    }
    size_t sum = 0;
    int ret;
    for (size_t i = 0; i < m; ++i) {
        ret = write_row_csv(f, n, mat[i]);
        if (ret < 0) {
            return ret;
        } else {
            sum += ret;
        }
    }
    ret = assert_size_csv(f, m, n);
    if (ret < 0) {
        return ret;
    }
    fclose(f);
    return sum;
}
