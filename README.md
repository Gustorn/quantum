QSpice
======

A quantum circuit simulator inspired by [simuleios's QSPICE](https://github.com/leios/QSPICE) library but using a completely different way of encoding qubits and quantum gates. It's written in Julia, a high performance programming language for numerical applications and provides a user friendly, text-based format for describing quantum circuits.

Features
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

Examples
--------

The following examples are provided with the library:

- Quantum Fourier Transform
- Shor's algorithm
- Superdense coding
- Quantum Teleportation
