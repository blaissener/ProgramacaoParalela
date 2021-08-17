#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

enum { SUCCESS, BAD_ARGUMENT };

// Function to sort a vector using bubble sort.
void bubble_sort(int *v, size_t n);

// Custom function to read command line arguments.
void read_arguments_or_abort(int argc, char *argv[], size_t *N, size_t *M);

/**
 * Read number of elements (N) from the command line, generate a
 * random vector of this size and sort it.  Repeat M times for
 * statistics (M is the second command line argument).
 */
int main(int argc, char *argv[]) {
  // Read N and M from command line.
  size_t N;
  size_t M;
  read_arguments_or_abort(argc, argv, &N, &M);

  // Initialize random number generator.
  srand(time(NULL));

  int *rnd_array = (int *)malloc(N * sizeof(int));

  // Time sorting M random arrays of N elements.
  double elapsed = 0;

  struct timeval t1, t2;

  for (size_t j = 0; j < M; j++) {

    // Generate new random values.
    for (size_t i = 0; i < N; ++i) {
      rnd_array[i] = ((double)rand() / RAND_MAX) * 1000.;
    }

    gettimeofday(&t1, NULL);

    bubble_sort(rnd_array, N);

    gettimeofday(&t2, NULL);

    double dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
    elapsed += dt;
  }

  // Show timing resuts.
  printf("%8zu %18.12f\n", N, elapsed / M);

  free(rnd_array);

  return SUCCESS;
}

/*
 * Implementation of the used functions.
 */

void bubble_sort(int *v, size_t n) {
  bool swapped;

  do {
    // No swapped elements yet.
    swapped = false;

    for (size_t i = 0; i < n - 1; ++i) {
      if (v[i + 1] < v[i]) {
        // Elements out of order. Swap!
        int tmp = v[i];
        v[i] = v[i + 1];
        v[i + 1] = tmp;

        // Register the swap.
        swapped = true;
      }
    }
    // One more element at the right position at the end.
    --n;
    // Repeat until no swaps (all elements in order).
  } while (swapped);

}

void read_arguments_or_abort(int argc, char *argv[], size_t *N, size_t *M) {
  // We need exactly three arguments (considering program name).
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <number of elements> <number of arrays>\n",
            argv[0]);
    exit(BAD_ARGUMENT);
  }

  // Read arguments.

  int n_read_OK = sscanf(argv[1], "%zu", N);
  if (n_read_OK != 1 || *N < 1) {
    fprintf(stderr, "Wrong argument for number of elements.\n");
    exit(BAD_ARGUMENT);
  }

  n_read_OK = sscanf(argv[2], "%zu", M);
  if (n_read_OK != 1 || *M < 1) {
    fprintf(stderr, "Wrong argument for number of repetitions.\n");
    exit(BAD_ARGUMENT);
  }
}
