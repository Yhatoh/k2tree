# Depth-First Representation of a $k^2$-Tree

This project implements a depth-first representation of a $k^2$-tree, an efficient data structure for representing sparse matrices. It leverages two external libraries: **sdsl-lite** for succinct data structures and **libsais** for efficient suffix array construction.

## Requirements

- **C++17** compiler
- **sdsl-lite** – choose one of the following versions:
  - [simongog's version](https://github.com/simongog/sdsl-lite)
  - [vgteam's version](https://github.com/vgteam/sdsl-lite) (includes support for rle-vector and more updated features)
- **libsais** – available at [libsais repository](https://github.com/IlyaGrebnov/libsais)

## Compilation

You can compile the project with or without optimization flags. By default, the Makefile assumes:
- **sdsl-lite** is located at `~/sdsl-lite`
- **libsais** is located in the `libsais` directory

### Default Build

```bash
make
```

### Release Build

To compile with optimization and release flags (`-O3 -m64 -DNDEBUG`):

```bash
make release
```

### Custom Library Paths

If you have custom paths for the libraries, you can override the defaults:

```bash
make SDSL_DIR=/your/path/to/sdsl-lite LIBSAIS_DIR=/your/path/to/libsais
```

## Usage

### Building the $k^2$-Tree

1. **Prepare Your Sparse Matrix**

   Create a text file representing your sparse matrix. The file should have the following format for a matrix of size $n \times n$ with $m$ ones:

   ```
   x1 y1
   x2 y2
   ...
   xm ym
   ```

2. **Build the $k^2$-Tree**

   Execute the build tool with the required parameters:

   ```bash
   path_to_folder/k2bp_build.x path_to_matrix/matrix.txt n m
   ```

   This command will generate a file named `matrix.txt.k2bp` in the same folder as your input matrix.

### Compressing the $k^2$-Tree

After building the $k^2$-tree, compress it using:

```bash
path_to_folder/k2bp_compr.x path_to_k2bp/matrix.txt.k2bp
```

This will produce a compressed tree file named `matrix.txt.k2bpi`.

### Operations

#### Information

To display information about the $k^2$-tree (e.g., space usage and tree size), run:

```bash
path_to_folder/k2bp_info.x path_to_k2bp/matrix.txt.k2bp
```

For the compressed $k^2$-tree (to show details about pruned interchangeable subtrees):

```bash
path_to_folder/k2bp_compr_info.x path_to_k2bp/matrix.txt.k2bpi
```

#### Multiplication

You can multiply two $k^2$-trees using the following commands:

- **For Uncompressed Trees:**

  ```bash
  path_to_folder/k2bp_mult_1.x path_to_k2bp1/matrix1.txt.k2bp path_to_k2bp2/matrix2.txt.k2bp
  ```

- **For Compressed Trees:**

  ```bash
  path_to_folder/k2bp_compr_mult.x path_to_k2bpi1/matrix1.txt.k2bpi path_to_k2bpi2/matrix2.txt.k2bpi
  ```

## Acknowledgments

- [sdsl-lite](https://github.com/simongog/sdsl-lite) and [vgteam's sdsl-lite](https://github.com/vgteam/sdsl-lite)
- [libsais](https://github.com/IlyaGrebnov/libsais)
