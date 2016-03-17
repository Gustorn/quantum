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

Netlist syntax:
===============

The library currently supports the following syntax for defining quantum circuits (see `test/test.qnl` for a full example):

```text
index:gate[connections](arguments)
```

|Syntax Element|Possible Values                                 |
|--------------|------------------------------------------------|
| `index`      | 64-bit signed integer                          |
| `gate`       | one of the possible gates<sup>[1]</sup>        |
| `connections`| comma-delimited `index` values                 |
| `arguments`  | comma-delimited arguments values<sup>[2]</sup> |

Where `connections` and `arguments` are optional, if the specific gate doesn't need it.

`[1]`: Gates and their arguments:

|Gate          | Arg 1                       | Arg 2                       | Arg 3                  |
|--------------|-----------------------------|-----------------------------|------------------------|
| `hadamard`   | `bit:      ` 1-based index  |                             |                        |
| `not`        | `bit:      ` 1-based index  |                             |                        |
| `cnot`       | `control:  ` 1-based index  | `flip:      ` 1-based index |                        |
| `ccnot`      | `control 1:` 1-based index  | `control 2: ` 1-based index | `flip: ` 1-based index |
| `swap`       | `bit 1:    ` 1-based index  | `bit 2:     ` 1-based index |                        |
| `cswap`      | `control:  ` 1-based index  | `bit 1:     ` 1-based index | `bit 2:` 1-based index |
| `sqrt_swap`  | `bit 1:    ` 1-based index  | `bit 2:     ` 1-based index |                        |
| `phase_shift`| `bit:      ` 1-based index  | `angle:     ` number        |                        |
| `pauli_x`    | `bit:      ` 1-based index  |                             |                        |
| `pauli_y`    | `bit:      ` 1-based index  |                             |                        |
| `pauli_z`    | `bit:      ` 1-based index  |                             |                        |
| `measure`    |                             |                             |                        |
| `partial_measure` | `bits:` any number of 1-based indices                |                        |
| `probe`      |                             |                             |                        |

Note: The quantum states are indexed in the following manner (using the `bra-ket` notation, replacing the individual qubits with their respective indices): `|1 2 3 4 5 6 ... >`

`[2]`: The possible arguments include:

* Integers and floating-point numbers
* Simple expressions containing at most one of the following operators: `+`, `-`, `*`, `\` and at most one (per side) of the following constants: `pi`, `i`. Some examples:
  * Valid: `1 + i`, `pi * i`, `1.0 * i`, `pi + 5.0i`, `pi / 8`
  * Invalid: `1 + 2 + 3`, `ipi`, `10 * pi * i`
