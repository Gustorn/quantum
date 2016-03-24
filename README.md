QSpice
======

A quantum circuit simulator inspired by [simuleios's QSPICE](https://github.com/leios/QSPICE) library. It's written in Julia, a high performance programming language for numerical applications.

Features:
---------

* Netlist parsing with a custom syntax
* The following quantum gates are implemented:
  * Hadamard
  * NOT: NOT, CNOT, CCNOT
  * SWAP: SWAP, CSWAP, SQRT_SWAP
  * Phase-Shift
  * Pauli-X/Y/Z
  * Measurement
  * Partial Measurement: both 1 and n-bit variants, in the computational basis
  * Apart from the standard gates, the library also includes a special "probe" gate that provides information about the quantum state

