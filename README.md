# Sparse Matrix Calculator in C

An interactive, menu-driven sparse-matrix calculator written in C99 as a
data-structures learning project. Only non-zero elements are stored, in a
linked structure that can be walked by row or by column. Up to ten matrices
(named `A`–`J`) can be held in memory at once.

This is a single-file educational implementation with a small regression test
suite. It is not a numerical library: see [Limitations](#limitations) before
relying on it for anything beyond small examples.

## Features

Available from the interactive menus:

- Create a matrix with given dimensions and non-zero entries, resize it, insert
  or update an element, delete an element, clear all of its elements.
- Add, subtract and multiply matrices; multiply by a scalar; transpose.
- Determinant and inverse (recursive cofactor expansion, so small matrices only).
- Full view (zeros shown), sparse view (non-zero entries with coordinates),
  dimensions only, and a list of all ten matrix slots.
- After a binary operation or an inverse, the result is printed as `R` and you
  are offered the chance to save it into any slot `A`–`J`.

Internal functions that are **not** exposed through the menus: `search` and
`copyMatrix` (copying is only reachable through the "save result" prompt).

## Repository layout

```text
sparse_matrix_calculator.c   the whole implementation, including main()
tests/test_regressions.c     regression tests (includes the .c file directly)
Makefile                     build, test and sanitizer targets
```

## Build and run

Requires a C99 compiler and the math library. Tested with GCC 13.3 on Ubuntu.

```bash
make
./matrix_calculator
```

Without Make:

```bash
gcc -std=c99 -Wall -Wextra -Wpedantic -o matrix_calculator sparse_matrix_calculator.c -lm
```

With GCC 13.3, both the calculator and the test driver compile with no
warnings under these flags.

## Usage

The main menu has three areas: **Matrix Management** (create, resize,
insert/update, delete, clear), **Matrix Operations** (command mode) and
**Display Options**.

### Conventions

- Matrix names are single **uppercase** letters `A`–`J`.
- Row and column indexes are **zero-based**.
- Matrix dimensions must be non-negative; negative sizes are rejected when
  creating or resizing.
- Elements are `float` (`matrix_entry`). Printed values use two decimals.
- Inserting at an existing coordinate updates the value.
- Inserting `0` never creates a stored element and does **not** overwrite an
  existing one. Use Delete Element to remove an entry.
- `scalar A 0` removes every stored element. An element whose value becomes
  exactly `0` after scalar multiplication (including by `float` underflow) is
  also removed.

### Input handling

Every prompt reads one whole line. Text that is not a valid number (or has
trailing characters) prints `Invalid numeric input.` and the prompt is shown
again, and leftover input no longer leaks into the next prompt or into the
first command. End of input exits the program cleanly.

### Entering elements when creating a matrix

The prompt reads all three values from one line. To store `1` at (0, 0), type:

```text
> 0 0 1
```

and type `-1` on a line by itself to finish. Entering `0 0 1` on one line is
accepted. (The Insert/Update Element menu item uses separate prompts for row,
column and value and is not affected.)

### Command mode

```text
add A B
subtract A B
multiply A B
transpose A
determinant A
inverse A
scalar A 2.5
exit
```

- `add`, `subtract`, `multiply` and `inverse` print the result as matrix `R`
  and ask whether to save it. They do not modify their operands.
- `transpose A` and `scalar A <value>` **modify `A` in place**.
- `determinant A` prints the value and changes nothing.
- Addition and subtraction require equal dimensions; multiplication requires
  `cols(left) == rows(right)`. This is checked even when an operand has no
  stored elements; an empty-operand product has the correct
  `rows(left) × cols(right)` size and no stored elements.

### Example

Create `A` as a 2×2 matrix with entries 1, 2, 3, 4 (Matrix Management →
Create Matrix, using the single-line entry format above), then in command mode:

```text
> Operation: determinant A
Determinant of A = -2.00
> Operation: inverse A
Matrix R [2 x 2]
 -2.00   1.00
  1.50  -0.50
```

## How it works

Each non-zero element is one node that sits in two sorted singly linked lists
at once: its row's list (via `right`) and its column's list (via `down`).
Row headers and column headers are themselves sorted linked lists, and a header
exists only while its row or column has at least one element.

```text
Non-zero elements: (0,1)  (0,4)  (2,1)

row headers:  [row 0] -> (0,1) -> (0,4)     via right pointers
              [row 2] -> (2,1)

col headers:  [col 1] -> (0,1) -> (2,1)     via down pointers
              [col 4] -> (0,4)
```

```c
typedef struct Sm_Node_Tag {
    matrix_entry data;                 /* float */
    int row, col;
    struct Sm_Node_Tag *right, *down;  /* next in row, next in column */
} Sm_Node;
```

Design notes that follow from the code:

- A node is shared between its row and column chains, so updating a value
  touches one node, and deleting an element unlinks it from both chains.
- Transpose reuses the existing element nodes: it swaps each node's
  `row`/`col` and `right`/`down`, and rebuilds only the header lists.
- Ten matrices live in a global registry (`NamedMatrix registry[10]`).

### Memory

On a 64-bit build, `sizeof(Sm_Node)` is 32 bytes, and each row or column header
is 24 bytes. A dense `float` matrix costs 4 bytes per cell, so ignoring headers
the sparse form only saves memory when fewer than roughly 1 in 8 cells is
non-zero. This is a derived estimate, not a benchmark.

### Time cost

There is no hashing or random access; everything walks linked lists.

- Insert, delete and search scan the row-header list, the column-header list,
  and the target row/column chains.
- Add, subtract and multiply traverse both operands, but build the result
  through `insertElement`, so they cost more than linear in the number of
  non-zeros.
- Determinant and inverse grow **factorially** with matrix size. Measured on
  one machine with an unoptimized build (`-O0`), dense determinant took about
  0.02 s at 8×8, 0.19 s at 9×9 and 2.0 s at 10×10; dense inverse took about
  0.19 s at 8×8. Sparse inputs are much cheaper because zero entries are skipped
  (a 12×12 diagonal matrix is effectively instant). Treat these as rough
  indications, not benchmarks.

## Tests

`tests/test_regressions.c` is deterministic and uses `assert`. It currently covers:

- Shrinking a matrix with `resizeMatrix` keeps the right elements, drops the
  rest, and leaves the row lists and column lists holding the same number of
  elements; `clearMatrix` empties both header lists.
- Fractional addition (`1.3 + 0.3`) and multiplication (`1.3 × 2.5`) to within 0.001.
- Updating an existing element succeeds and changes the stored value.
- Inverting a 1×1 matrix.
- Negative and out-of-range coordinates are rejected by `insertElement`, and a
  negative size passed to `resizeMatrix` is ignored.
- Inserting `0` does not overwrite an existing value; scalar multiplication by
  `0` empties both indexes.
- Addition and multiplication reject mismatched dimensions even when an operand
  is empty, and an empty-operand product has the expected dimensions.
- Small diagonal matrices and mixed-scale diagonal matrices retain their
  non-zero determinants and inverse cofactors.

Not covered: delete, search, subtract, transpose, determinant or inverse beyond
1×1, copy, command parsing, the input-reading helpers and interactive menus,
and allocation failure.

```bash
make matrix_regressions
./matrix_regressions
```

### Sanitizers

Build and run the same tests with AddressSanitizer and UndefinedBehaviorSanitizer:

```bash
make asan
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1 ./matrix_regressions_asan
```

With GCC 13.3 on Ubuntu this completes with no sanitizer findings. That result
applies to what the tests exercise, not to the whole program. On Windows, run
these commands under WSL2.

(The out-of-bounds message printed several times during the run comes from the
invalid-coordinate tests and is expected.)

## Changes from the original version

Compared with the first version of this repository:

- Addition and multiplication accumulated results in an `int`, truncating
  fractional values; they now use `float`.
- Updating an existing element printed "Unexpected duplication" and returned
  failure; it now succeeds.
- `clearMatrix` and `resizeMatrix` read or freed nodes after they had been
  deleted, and left column chains pointing at freed memory; they now save the
  next pointer first and go through `deleteElement`.
- Negative coordinates and sizes are rejected, node allocation is checked, and
  1×1 inverse has an explicit case.
- Add and multiply validate dimensions up front, including for empty operands.
- Scalar multiplication removes elements that become zero instead of storing
  zero-valued nodes.
- Prompts read whole lines and validate numbers instead of using bare `scanf`.
- Singularity is judged with a tolerance scaled to the product of the matrix
  row magnitudes, and small non-zero cofactors are retained.
- The build is now `-std=c99 -Wall -Wextra -Wpedantic` with a Makefile, tests
  and sanitizer target, and compiles without warnings on GCC 13.3.

## Limitations

Behaviors confirmed by running the program:

- End of input exits the program cleanly.
- Matrix creation accepts the documented single-line `row col value` format.
- Singularity detection is scaled to the matrix magnitude, and non-zero
  cofactors are not discarded merely because they are small. Results remain
  subject to `float` precision and recursive cofactor-expansion limits.
- Determinant and inverse become impractical beyond roughly 10×10 for dense
  matrices (see [Time cost](#time-cost)). Results are subject to `float`
  rounding.

From reading the code:

- Allocation failures are handled in the node-creation and insert path, but not
  uniformly across every operation, and that handling is untested.
- There is a global registry of ten matrices, uppercase names only, and no file
  input or output.