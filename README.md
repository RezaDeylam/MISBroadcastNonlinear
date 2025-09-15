# MISBroadcastNonlinear
MIS in Linear Function Computation Broadcast


This project contains MATLAB code to construct **confusion graphs** for a 3-user function computation broadcast model (q=3) and to enumerate all **maximal independent sets (MIS)**.

Although the original example is linear, this repository is structured so that it can be extended to **nonlinear broadcast functions** as well.

## Features
- Build the **union confusion graph** from broadcast functions.
- Correct adjacency rule: two states are adjacent if `Xi` is equal **and** `fi` differs.
- Enumerate all **maximal independent sets (MIS)** using the Bron–Kerbosch algorithm with pivoting.
- Verified on the 3-user linear example: the union graph has **225 MIS**.
- Easily adaptable to nonlinear functions.

## Requirements
- MATLAB R2021a or later (most versions will work).
- No extra toolboxes required.

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/MISBroadcastNonlinear.git
   cd MISBroadcastNonlinear/src
