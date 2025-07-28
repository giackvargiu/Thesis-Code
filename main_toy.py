#!/usr/bin/env python3

import random

def generate_grid(rows=10, cols=10, min_val=0, max_val=9):
    """
    Generate a 2D grid of random integers.

    Parameters:
    - rows (int): Number of rows in the grid.
    - cols (int): Number of columns in the grid.
    - min_val (int): Minimum random integer (inclusive).
    - max_val (int): Maximum random integer (inclusive).

    Returns:
    - list of lists containing the generated grid.
    """
    grid = []
    for i in range(rows):
        row = []
        for j in range(cols):
            val = random.randint(min_val, max_val)
            row.append(val)
        grid.append(row)
    return grid


def print_grid(grid):
    """
    Print the 2D grid to the terminal, row by row.
    """
    for row in grid:
        print(' '.join(str(val) for val in row))


def main():
    # Generate a 10x10 grid of random integers between 0 and 9
    grid = generate_grid()
    print("Generated 10x10 grid of random integers (0-9):")
    print_grid(grid)


if __name__ == "__main__":
    main()
