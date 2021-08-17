#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

enum { SUCCESS, BAD_ARGUMENT };

// Hand-mplemented quicksort of a vector of ints (for comparison).
void quick_sort(int *v, size_t n);

// Custom function to read command line arguments
void read_arguments_or_abort(int argc, char *argv[], size_t *N, size_t *M);

/**
 * Read number of elements (N) from the command line, generate a
 * random vector of this size and sort it.
 * Repeat M times for statistics (M is the second command line argument).
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

    quick_sort(rnd_array, N);

    gettimeofday(&t2, NULL);

    double dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
    elapsed += dt;
  }

  // Show timing resuts.
  printf("%8zu %18.12f\n", N, elapsed / M);

  free(rnd_array);

  return SUCCESS;
}

/**
 * Auxiliary recursive function.
 *
 * Sorts recursively elements from ini (included) to fin (excluded) using
 * quicksort.
 */
void quick_sort_rec(int *ini, int *fin) {
  // See if there is something to sort.
  ptrdiff_t interval_size = fin - ini;
  if (interval_size > 2) {
    // More than two elements: needs sorting.

    // Choose first element as key.
    int key = *ini;

    // Split elements:
    // 1) [ini, smaller_border) are elements less than or equal to key.
    // 2) [larger_border, fin) are elements greater than key.
    // 3) [smaller_border, larger_border) are yet unchecked elements.
    int *smaller_border = ini + 1;
    int *larger_border = fin;
    while (smaller_border != larger_border) {
      // There are yet unknowns.
      // Put value pointed to by smaller_border at the right place.
      if (*smaller_border <= key) {
        // It's already in the right place. Only adjust border.
        ++smaller_border;
      } else {
        // Swap *smaller_border with *(larger_border - 1), i.e. the last unknown.
        int tmp = *smaller_border;
        *smaller_border = *(larger_border - 1);
        *(larger_border - 1) = tmp;

        // One more large value in position.
        --larger_border;
      }
    }

    // Put the key element in the division point.
    *ini = *(smaller_border - 1);
    *(smaller_border - 1) = key;

    // Recursive call for the two subsequences.
    quick_sort_rec(ini, smaller_border - 1);
    quick_sort_rec(larger_border, fin);
  } else if (interval_size == 2 && *ini > *(ini + 1)) {
    // The two elements are out of order. Swap.
    int tmp = *ini;
    *ini = *(ini + 1);
    *(ini + 1) = tmp;
  }
}


/**
 * Sort all elements of array v using quick sort.
 */
void quick_sort(int *v, size_t n) { quick_sort_rec(v, v + n); }


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
