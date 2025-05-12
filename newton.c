#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <getopt.h>
#include <threads.h>
#include <stdatomic.h>

#define PI acos(-1.0)
#define COLOR_STR_SIZE 15

typedef struct {
    int lval;
    int power;
    double step;
    double complex *roots;
    int *root_img;
    int *iter_img;
    int start_row;
    int end_row;
    int max_iter;
    atomic_int *row_ready;
} thrd_data;

typedef struct {
    int lval;
    int power;
    int *root_img;
    int *iter_img;
    atomic_int *row_ready;
    char **attractor_colors;
    char **convergence_colors;
    char *divergence_color;
} writer_data;

void parse_args(int argc, char *argv[], int *tval, int *lval, int *power) {
    *tval = 1;     
    *lval = 1000;  
    *power = -1;   

    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "-t", 2) == 0){
            *tval = atoi(argv[i] + 2);
        } else if (strncmp(argv[i], "-l", 2) == 0){
            *lval = atoi(argv[i] + 2);
        } else if (argv[i][0] != '-') {
            *power = atoi(argv[i]);
        } else {
            fprintf(stderr, "unknow arg: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if (*power <= 0) {
        fprintf(stderr, "error: degree must be int > 0.\n");
        exit(EXIT_FAILURE);
    }
}

void unity_roots(int d, double complex roots[]) {
    double theta;
    for (int k = 0; k < d; k++) {
        theta = 2 * PI * k / d;
        roots[k] = cos(theta) + sin(theta) * I;
    }
}

double complex complex_pow(double complex base, int exponent) {
    double complex result = 1.0 + 0.0 * I;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    return result;
}

void newton(int power, const int max_iter,double complex roots[], double complex *z, int *n_iter, int *conv_root_idx) {
    double tol = 1e-3;
    double tol2 = tol * tol;
    int i;
    double complex z_prev;
    for (i = 0; i < max_iter; i++) {
        z_prev = *z;
        // Check for divergence
        double z_abs2 = creal(*z) * creal(*z) + cimag(*z) * cimag(*z);
        if (z_abs2 > 1e20) {
            *n_iter = -1;
            *conv_root_idx = -1;
            return;
        }
        // Compute z^{n-1} using exponentiation by squaring
        double complex z_n_minus1 = complex_pow(*z, power - 1);
        double complex z_n = z_n_minus1 * (*z);
        // Compute f(z) and f'(z)
        double complex f_z = z_n - 1;
        double complex df_z = power * z_n_minus1;
        // Newton iteration step
        double a_re = creal(f_z);
        double a_im = cimag(f_z);
        double b_re = creal(df_z);
        double b_im = cimag(df_z);
        double denom = b_re * b_re + b_im * b_im;
        double real_part = (a_re * b_re + a_im * b_im) / denom;
        double imag_part = (a_im * b_re - a_re * b_im) / denom;
        *z -= real_part + I * imag_part;
        // Check for convergence
        double delta_re = creal(*z) - creal(z_prev);
        double delta_im = cimag(*z) - cimag(z_prev);
        if (delta_re * delta_re + delta_im * delta_im < tol2) {
            *n_iter = i;
            break;
        }
    }
    if (i == max_iter) {
        *n_iter = -1;
        *conv_root_idx = -1;
        return;
    }
    // Determine which root z converged to
    double min_dist2 = 1e20;
    int closest_root = -1;
    for (int k = 0; k < power; k++) {
        double delta_re = creal(*z) - creal(roots[k]);
        double delta_im = cimag(*z) - cimag(roots[k]);
        double dist2 = delta_re * delta_re + delta_im * delta_im;
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closest_root = k;
        }
    }
    if (min_dist2 < tol2) {
        *conv_root_idx = closest_root;
    } else {
        *conv_root_idx = -1; // converged but not to a root
    }
}


int color_mapper_root(int root_idx, int n_roots) {
    int step = 255 / n_roots;
    int r = ((root_idx) * step) % 256;
    int g = ((root_idx + 1) * step) % 256;
    int b = ((root_idx + 2) * step) % 256;
    return (r << 16) + (g << 8) + b;
}

int color_mapper_iter(int n_iter) {
    int r = (n_iter % 256);
    int g = (n_iter * 5 % 256);
    int b = (n_iter * 10 % 256);
    return (r << 16) + (g << 8) + b;
}

int thread_ctrl(void *arg) {
    thrd_data *data = (thrd_data *)arg;
    for (int i = data->start_row; i < data->end_row; i++) {
        for (int j = 0; j < data->lval; j++) {

            double Re = -2.0 + j * data->step;
            double Im = -2.0 + i * data->step;
            double complex z = Re + Im * I;
            int n_iter, conv_root_idx;
            newton(data->power, data->max_iter, data->roots, &z, &n_iter, &conv_root_idx);

            int idx = i * data->lval + j;
            if (n_iter != -1) {
                if (conv_root_idx != -1) {
                    data->root_img[idx] = conv_root_idx; // storing root idx 
                } else {
                    data->root_img[idx] = -1; // divergence
                }
                data->iter_img[idx] = n_iter;
            } else {
                data->root_img[idx] = -1; // divergence 
                data->iter_img[idx] = 0;
            }
        }
        atomic_store(&data->row_ready[i], 1); // storing which row is ready
    }
    return 0;
}

int write_thread(void *arg) {
    writer_data *wdata = (writer_data *)arg;
    int lval = wdata->lval;
    int *root_img = wdata->root_img;
    int *iter_img = wdata->iter_img;
    atomic_int *row_ready = wdata->row_ready;
    char **attractor_colors = wdata->attractor_colors;
    char **convergence_colors = wdata->convergence_colors;
    char *divergence_color = wdata->divergence_color;

    // Open files
    char attractor_filename[256];
    char convergence_filename[256];
    sprintf(attractor_filename, "newton_attractors_x%d.ppm", wdata->power);
    sprintf(convergence_filename, "newton_convergence_x%d.ppm", wdata->power);

    FILE *attr_file = fopen(attractor_filename, "w");
    FILE *conv_file = fopen(convergence_filename, "w");

    if (!attr_file || !conv_file) {
        fprintf(stderr, "Error opening output files.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(attr_file, "P3\n%d %d\n255\n", lval, lval);
    fprintf(conv_file, "P3\n%d %d\n255\n", lval, lval);

    size_t buffer_size = lval * COLOR_STR_SIZE;
    char *attr_row_buffer = malloc(buffer_size);
    char *conv_row_buffer = malloc(buffer_size);

    for (int i = 0; i < lval; i++) {
        while (atomic_load(&row_ready[i]) == 0) {
            thrd_yield();
        }
        size_t attr_offset = 0;
        size_t conv_offset = 0;

        for (int j = 0; j < lval; j++) {
            int idx = i * lval + j;

            int attractor_index = root_img[idx];
            const char *color_str;
            if (attractor_index >= 0 && attractor_index < wdata->power) {
                color_str = attractor_colors[attractor_index];
            } else {
                color_str = divergence_color;
            }
            strcpy(&attr_row_buffer[attr_offset], color_str);
            attr_offset += strlen(color_str);

            int n_iter = iter_img[idx];
            if (n_iter < 0 || n_iter > 50) {
                n_iter = 50;
            }
            color_str = convergence_colors[n_iter];
            strcpy(&conv_row_buffer[conv_offset], color_str);
            conv_offset += strlen(color_str);
        }

        fwrite(attr_row_buffer, sizeof(char), attr_offset, attr_file);
        fwrite("\n", sizeof(char), 1, attr_file);
        fwrite(conv_row_buffer, sizeof(char), conv_offset, conv_file);
        fwrite("\n", sizeof(char), 1, conv_file);
    }
    free(attr_row_buffer);
    free(conv_row_buffer);
    fclose(attr_file);
    fclose(conv_file);
    return 0;
}

int main(int argc, char **argv) {
    int tval, lval, power;

    parse_args(argc, argv, &tval, &lval, &power);

    if (power <= 0) {
        fprintf(stderr, "Error: power must be a positive integer.\n");
        return EXIT_FAILURE;
    }

    // allocate array for roots
    double complex *roots = malloc(sizeof(double complex) * power);
    unity_roots(power, roots);

    // allocate array for images
    int *root_img = calloc(lval * lval, sizeof(int));
    int *iter_img = calloc(lval * lval, sizeof(int));

    // array to keep track of processed rows
    atomic_int *row_ready = malloc(sizeof(atomic_int) * lval);
    for (int i = 0; i < lval; i++) {
        atomic_init(&row_ready[i], 0);
    }

    double step = 4.0 / (lval - 1);
    const int max_iter = 10000;

    char **attractor_colors = malloc(sizeof(char *) * power);
    for (int i = 0; i < power; i++) {
        int color = color_mapper_root(i, power);
        int r = (color >> 16) & 0xFF;
        int g = (color >> 8) & 0xFF;
        int b = color & 0xFF;
        attractor_colors[i] = malloc(COLOR_STR_SIZE); // 12 bytes for color plus NULL terminator
        sprintf(attractor_colors[i], "%d %d %d ", r, g, b);
    }

    char divergence_color[12];
    sprintf(divergence_color, "0 0 0 ");

    int max_iterations = 50;
    char **convergence_colors = malloc(sizeof(char *) * (max_iterations + 1));
    for (int i = 0; i <= max_iterations; i++) {
        int color = color_mapper_iter(i);
        int r = (color >> 16) & 0xFF;
        int g = (color >> 8) & 0xFF;
        int b = color & 0xFF;
        convergence_colors[i] = malloc(COLOR_STR_SIZE);
        sprintf(convergence_colors[i], "%d %d %d ", r, g, b);
    }

    writer_data wdata;
    wdata.lval = lval;
    wdata.power = power;
    wdata.root_img = root_img;
    wdata.iter_img = iter_img;
    wdata.row_ready = row_ready;
    wdata.attractor_colors = attractor_colors;
    wdata.convergence_colors = convergence_colors;
    wdata.divergence_color = divergence_color;

    thrd_t writer_thread;
    thrd_create(&writer_thread, write_thread, &wdata);

    thrd_t *compute_threads = malloc(sizeof(thrd_t) * tval);
    thrd_data *thread_data = malloc(sizeof(thrd_data) * tval);

    int rows_per_thread = lval / tval;
    int extra_rows = lval % tval;
    int row_start = 0;

    for (int i = 0; i < tval; i++) {
        int rows = rows_per_thread + (i < extra_rows ? 1 : 0);
        thread_data[i].lval = lval;
        thread_data[i].power = power;
        thread_data[i].step = step;
        thread_data[i].roots = roots;
        thread_data[i].root_img = root_img;
        thread_data[i].iter_img = iter_img;
        thread_data[i].start_row = row_start;
        thread_data[i].end_row = row_start + rows;
        thread_data[i].max_iter = max_iter;
        thread_data[i].row_ready = row_ready;
        row_start += rows;
        thrd_create(&compute_threads[i], thread_ctrl, &thread_data[i]);
    }
    // wait for compute threads to finish
    for (int i = 0; i < tval; i++) {
        thrd_join(compute_threads[i], NULL);
    }
    // wait for writer thread to finish
    thrd_join(writer_thread, NULL);
    // clean up allocated memory
    for (int i = 0; i < power; i++) {
        free(attractor_colors[i]);
    }
    free(attractor_colors);

    for (int i = 0; i <= max_iterations; i++) {
        free(convergence_colors[i]);
    }
    free(convergence_colors);
    free(root_img);
    free(iter_img);
    free(row_ready);
    free(compute_threads);
    free(thread_data);
    free(roots);

    return 0;
}


