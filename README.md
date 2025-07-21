# PhaseTransitions

A Julia package for simulating and analyzing percolation and phase transitions on lattices.

## Features
- Lattice generation
- Percolation simulation
- Visualization tools
- Cluster analysis

## Installation
Clone the repository and use Julia's package manager:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Usage Example
Here is a basic example of how to use the package:

```julia
using PhaseTransitions

# Example: Run a percolation simulation
lattice = PhaseTransitions.generate_lattice(50, 50)
perc = PhaseTransitions.site_percolation(lattice, p=0.6)
PhaseTransitions.visualize(perc)
```

See the `src/` and `test/` directories for more examples and details.

## Contributing
Pull requests and issues are welcome!

## License
MIT License 