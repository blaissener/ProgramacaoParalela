#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum Error { SUCCESS, BAD_ARGUMENT };

double dot_product(const double *a, const double *b, size_t N);

double random_uniform() { return (double)rand() / (RAND_MAX + 1.); }

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Give the size of the vectors and the number of threads in "
                    "the command line.\n");
    return BAD_ARGUMENT;
  }

  size_t N;
  int n_read_OK = sscanf(argv[1], "%zu", &N);

  if (n_read_OK != 1 || N < 1) {
    fprintf(stderr, "The number of elements must be positive.\n");
    return BAD_ARGUMENT;
  }

  int num_threads;
  n_read_OK = sscanf(argv[2], "%d", &num_threads);

  if (n_read_OK != 1 || num_threads <= 0) {
    fprintf(stderr, "The number of threads must be positive.\n");
    return BAD_ARGUMENT;
  }

  // Random number generation
  srand(time(NULL));

  double *a = (double *)malloc(N * sizeof(double));
  double *b = (double *)malloc(N * sizeof(double));

  for (size_t i = 0; i < N; ++i) {
    a[i] = random_uniform();
    b[i] = random_uniform();
  }

  omp_set_num_threads(num_threads);

  double dotprod = dot_product(a, b, N);

  printf("The dot product is %f.\n", dotprod);

  return SUCCESS;
}

double dot_product(const double *a, const double *b, size_t N) {
  double s = 0.0;

#pragma omp parallel for default(none) shared(a, b, N) reduction(+:s)     \
    schedule(static, 10)
  for (size_t i = 0; i < N; ++i) {
    s += a[i] * b[i];
  }

  return s;
}
