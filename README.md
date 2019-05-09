# The Xyce&trade; Parallel Electronic Simulator

[![Xyce](doc/Common_Guide_Files/xyce_flat_white.png)](https://xyce.sandia.gov)

## About Xyce

[Xyce](https://xyce.sandia.gov) (z&#x012B;s, rhymes with "spice") is an open
source, SPICE-compatible, high-performance analog circuit simulator, capable of
solving extremely large circuit problems by supporting large-scale parallel
computing platforms. It also supports serial execution on all common desktop
platforms, and small-scale parallel runs on Unix-like systems. In addition to
analog electronic simulation, Xyce has also been used to investigate more
general network systems, such as neural networks and power grids. In providing
an Open Source version of Xyce to the external community, Sandia hopes to
contribute a robust and modern electronic simulator to users and researchers in
the field.

While designed to be SPICE-compatible, Xyce is not a derivative of SPICE. Xyce
was designed from scratch to be a parallel simulation code, written in C++ and
using a message-passing implementation (MPI). Xyce also leverages Sandia's
open-source solver library, [Trilinos](https://github.com/trilinos/Trilinos),
which includes a number of circuit-specific solvers, such as the KLU direct
solver. With its modular and flexible design, Xyce applies abstract interfaces
to enable easy development of different analysis types, solvers and models.

Xyce is compatible with SPICE-based codes, in that it supports a canonical set
of SPICE compact models and standard SPICE analysis methods, such as
steady-state (`.DC`), transient (`.TRAN`), small-signal frequency domain
(`.AC`), and noise (`.NOISE`). However, Xyce goes beyond most SPICE-based codes
in a number of ways, including support for a large number of non-traditional
models, such as neuron and reaction network models. Xyce also supports Harmonic
Balance analysis (`.HB`), random sampling analysis, sensitivity calculations,
and post processing of the simulation metrics (`.FOUR` and `.MEASURE`).

## Binaries, building and more

For more information about Xyce&mdash;including documentation, build
instructions and binary executables&mdash;see the [Xyce
Homepage](https://xyce.sandia.gov) at [Sandia National
Laboratories](https://www.sandia.gov).

## Contributing

We welcome bug reports and enhancement requests. Those can be done through the
GitHub "Issues" area, or [other
methods](https://xyce.sandia.gov/contact_us.html). Due to internal
restrictions, however, it is difficult for us to accept external contributions
at this time, which includes pull requests. Nevertheless, if you would like to
discuss the possibility of a contribution, please contact us.

## Copyright and license

Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

Xyce&trade; is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Xyce&trade; is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

A copy of the GNU General Public License is included in the
[COPYING](./COPYING) file, or see <http://www.gnu.org/licenses/>.

## Acknowledgements

Xyce has been funded by the NNSA's Advanced Simulation and Computing (ASC)
Campaign, the DARPA POSH program, and the Laboratory Directed Research and
Development program at Sandia National Laboratories. Sandia National
Laboratories is a multimission laboratory managed and operated by National
Technology & Engineering Solutions of Sandia, LLC, a wholly owned subsidiary of
Honeywell International Inc., for the U.S. Department of Energy's National
Nuclear Security Administration under contract DE-NA0003525.

SAND2019-5105 O
