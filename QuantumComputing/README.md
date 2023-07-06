# Quantum Computing with Qiskit

<p align="center">
  <img src="https://miro.medium.com/v2/resize:fit:810/0*6qh01lw6kdaC73KF.png">
</p>

This repository is a collection of Jupyter notebooks aimed at giving a comprehensive understanding of Quantum Computing using Qiskit, an open-source quantum computing framework developed by IBM. It explores creating quantum circuits, applying different gates, advanced circuit techniques and much more, all using the beautiful language of quantum mechanics.

## Notebooks

1. **1_Circuit_Basics.ipynb** - Initiation to the quantum world, introducing quantum circuits. Here you'll start from the famous Dirac's bra-ket notation (|ψ⟩) and move to actual implementation of quantum circuits using Qiskit.

2. **2_1_Circuits_Operations.ipynb** - Detailed coverage of single-qubit operations. Here, you'll see how Pauli gates (X, Y, Z), Clifford gates, and more are represented in the quantum world, and how to use them in Qiskit. Each gate corresponds to a unitary transformation, e.g. X-gate is represented as:

```latex
    X = [0 1]
        [1 0]
```
3. **2_2_Circuits_Operations.ipynb** - Explores duo qubit operations which are fundamental for entanglement and other quantum properties. For example, the CNOT gate represented as:

```latex
    CNOT = [1 0 0 0]
           [0 1 0 0]
           [0 0 0 1]
           [0 0 1 0]
```
4. **2_3_Circuits_Operations.ipynb** - Discusses non-unitary operations like Measurements, Conditional Operations, and Initialization.

5. **3_1_Advanced_Circuit.ipynb** - Discusses advanced circuit techniques such as opaque gates, composite gates, and parameterized circuits.

6. **3_2_Advanced_Circuit.ipynb** - Dives deep into the various aspects of operators in quantum computing.

7. **4_Quantum_Fourier_Transform.ipynb** - Discusses and implements the Quantum Fourier Transform (QFT), the quantum analogue of the classical Fourier Transform, but exponentially faster!


## Requirements

- Python 3.7 or above
- Jupyter Notebook
- Qiskit

## Installation

To install the Qiskit package, run the following command:

```bash
pip install qiskit
```

Or you can install via conda:

```bash
conda install -c conda-forge qiskit
```

Launch Jupyter notebook:

```bash
jupyter notebook
```

## Contributing

Contributions are welcome.

## Contact

Please feel free to contact the maintainers of this repository for questions, suggestions or feedback. Quantum world is fascinating, let's explore it together!
