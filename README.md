# normal_corner

This code produces a corner plot for analytical multi-dimensional Gaussian distribution, using covariance matrix and mean matrix. It also allows us to plot another distribution, with reduced dimensionality, on top. I.e. a distribution for a case where we fixed one variable.

## Examples

It's as simple as `figure_1 = normal_corner.normal_corner(covariance_matrix,mean_vector,variable_labels)`!

![Analytical Gaussian corner plot](https://github.com/bvgoncharov/normal_corner/blob/master/example_1.png "Analytical Gaussian corner plot")

See demo.py for detailed usage examples.

## Installation

Option 1:
 - `pip install normal_corner`

Option 2:
 - `python setup.py install` in cloned directory

## Documentation

The main component is a normal\_corner function inside a normal\_corner package.

```python
def normal_corner(covm,mean,varlabels,fixedvarindex=None,fixedvarvalue=None,
           covm2=None,mean2=None,scale_factor=3,diagnostic=False,
           color='red',color2='blue', fig=None, **fig_kw):
```

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
 - `fig` : Matplotlib figure to plot to, possibly output of corner.corner to plot on top of MCMC corner plot
 - additional keywords are passed to figure().

