# normal_corner

This code produces a corner plot for analytical multi-dimensional Gaussian distribution, using covariance matrix and mean matrix. It also allows us to plot another distribution, with reduced dimensionality, on top. I.e. a distribution for a case where we fixed one variable.

See demo.py for examples.

## Installation

Option 1:
 - `pip install normal_corner`

Option 2:
 - `python setup.py install` in cloned directory

## Documentation

The main component is function normal\_corner inside a normal\_corner package.
Below is a description of inputs and outputs.

#### Output:
 - Matplotlib figure object with a corner plot

#### Main input:
 - `covm` : covariance matrix, _numpy array_, NxN.
 - `mean` : mean matrix, _numpy array_, 1xN.
 - `varlabels` : labels for plotting, 1xN, _list of str_ in LaTex format, between ($$).

#### Input for a second distribution on top:
 - `fixedvarindex` : index of variable that we do not use (fix), _int_, starting from 0. If not _None_, covm2 and mean2 must not be _None_.
 - `fixedvarvalue` : value of fixed variable, _float_. Leave None not to plot fixed value.
 - `newcov` and `newmean` : new covariance and mean matrices, same format, as above.

#### Optional input:
 - `scale_factor` : scale factor for plotting area, _float_, in sigma.
 - `diagnostic` : an option to print out some diagnostic messages, _bool_.
 - `color` : color for a main Normal distribution, _str_.
 - `color2` : color for a secondary Normal distribution, with reduced dimensionality, _str_.
