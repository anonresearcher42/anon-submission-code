# Anonymous Code for Influence Maximization Research

This repository provides the source code for an anonymous submission related to:

**RFIM**: multiple community structures.
(Manuscript currently under review; details are intentionally omitted.)

## Base Code Acknowledgement

This code builds upon the implementation of:

Qintian Guo, Sibo Wang, Zhewei Wei, and Ming Chen. 2020. Influence Maximization Revisited: Efficient Reverse Reachable Set Generation with Bound Tightened. In SIGMOD.

## Compilation

To compile the code, use the following command:

    g++ -std=c++17 -O3 rfim.cpp -o rfim

## Usage

### Step 1: Format the Graph

Before running the algorithm, format the graph:

    ./rfim -func=format -graphname=facebook -pdist=wc

### Step 2: Run the Algorithm

Execute the sg-hist method with the specified parameters:

    ./rfim -func=rfim -graphname=facebook -pdist=wc -seedsize=100 -eps=0.1 -delta=0.05 -method=sgh

## Important Parameter

Option       | Type    | Description
------------ | ------- | --------------------------------------------
func         | string  | Task to perform: `format` or `rfim`
graphname    | string  | Name of the input graph
pdist        | string  | Propagation model; we use `wc`
seedsize     | int     | Number of seeds to select
eps          | double  | Error bound for approximation
delta        | double  | failure probability
method       | string  | Method to use: `sgh`, `sg`, or `ag`

**Method options:**
- `sgh`: SG-HIST (proposed method)
- `sg`: SingleGreedy (baseline method)
- `ag`: AllGreedy (baseline method)

## Related Projects

For reference and comparison, please see the following open-source projects:

- SUBSIM / HIST: https://github.com/qtguo/subsim
- OPIM / OPIM-C: https://github.com/tangj90/OPIM