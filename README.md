# hsc
![example workflow](https://github.com/fraterenz/hsc/actions/workflows/clippy-fmt.yml/badge.svg)
![workflow](https://github.com/fraterenz/hsc/actions/workflows/test.yml/badge.svg)
[![codecov](https://codecov.io/gh/fraterenz/hsc/branch/master/graph/badge.svg?token=6JTE5AV0QO)](https://codecov.io/gh/fraterenz/hsc)

Stochastic simulations of the genetic changes in blood occurring during ageing.

## HSC dynamics during ageing

In the bone marrow, hematopoietic stem cells (HSCs) maintain lifelong blood production (hematopoiesis).  
In young, healthy individuals, many HSCs contribute roughly equally to this process, keeping the blood genetically diverse.

As we age, however, fewer HSCs actively contribute to hematopoiesis. Some HSCs acquire advantageous mutations that allow them to proliferate faster than others, leading to clonal competition and a reduction in genetic diversity in the blood.

`hsc` simulates this age-related loss of genetic heterogeneity assuming competing HSC clones. `hsc` simulates HSCs and mutations, and tracks over time:

- the HSCs in the system
- the neutral mutations carried by each HSC
- the positively selected mutations defining fit clones

See [here](output.md) for more information about the output of the simulations and [here](https://fraterenz.github.io/hsc/hsc/) for the documentation.

## Usage
`hsc` can be ran from the command line: `hsc [OPTIONS] <DIR> <COMMAND>`, see also `hsc --help` and `hsc -h`.
There are two commands `<COMMAND>`, which are `moran` and `exp-moran`.
They have their options as well, see `hsc <COMMAND> -h`.

See [here](output.md) for more information about the output of the simulations and [here](https://fraterenz.github.io/hsc/hsc/) for the documentation.

Use the python package [hscpy](https://github.com/fraterenz/hscpy) to plot and process the output of the simulations.
