# IBD Generator

A Rust library for generating Incomplete Block Designs (IBDs) with prohibitions. This is a port of the IBD generation code from the AlgDesign R package.

The algorithm creates an optimal design based on D-optimality. For smaller combinations of v, n_b, k, the algorithm should produce a BIBD (if one exists, and presuming no prohibitions of course).

## Features

- Generate optimal incomplete block designs (IBDs)
- Support for prohibited treatment pairs
- Configurable parameters:
  - Number of treatments (v)
  - Number of blocks (n_b) 
  - Block size (k)
  - Number of iterations
  - Number of repeats per iteration

## Library Usage

To use the library, you can call the `ibdgen` function with the desired parameters. Here's an example:

```rust
let result = ibdgen(v, n_b, block_size, iter, n_repeats, prohibited_pairs);
```
where:
- v is the number of treatments total
- n_b is the number of blocks
- block_size is the number of treatments per block
- n_repeats is the number of repeats per iteration within a single run. Default is 5.
- iter is the number of times the whole process is repeated. For smaller designs (e.g. 7,7,3). The default of 10 is plenty to find the best design. Otherwise - experiment.
- prohibited_pairs is a vector of tuples representing the prohibited pairs. For example, if you want to prohibit the pair (1,2), you would pass vec![("1,2".to_string())].

This will generate an IBD with the specified parameters and return the best solution found.

## Command Line Usage
The command line interface is a simple wrapper around the ibdgen function. 

Usage:
```
ibdgen 9 12 3 --n_repeats 5 --iter 10 --prohibited_pairs "1,2,3,4"
```
creates an IBD with v=9, n_b=12, k=3 design with prohibited pairs (1,2) and (3,4). --n_repeats and --iter are optional, and default to 5 and 10 respectively.

## License
This project is licensed under the GPL3 License - see the LICENSE file for details.

## Notes
The code works, but neither my Rust nor my linear algebra are world-class, so doubtless many improvements are possible. In particular, the code is slow for large combos of v, n_b, k given the find_delta_block function which iterates over all pairs of blocks several times.

## TODO
- Make the code faster
- Much better error handling. In particular, there is very little bounds checking.
- Needs some serious refactoring
- Better documentation
- Better tests
