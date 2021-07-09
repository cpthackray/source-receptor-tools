# source-receptor-tools

<br />

[(See included PFAS example)](pfas_example.ipynb)
### NOTE: Included data is for illustrative purposes only!
<br />
<br />

# Requirements
 - [xarray](http://xarray.pydata.org/en/stable/getting-started-guide/installing.html)
 - [numpy](https://numpy.org/install/)
 
<br />

# Metadata Structure

Based on tree structure:

`Source Type -> Chemical Species -> Exposure Pathway -> Spatiotemporal Information`

Relies on NetCDF/xarray for handling the spatiotemporal dimensions.

Implemented with \*.nc files using the Variable convention:
```python
name = f'{SOURCETYPE}_{CHEMICALSPECIES}_{EXPOSUREPATHWAY}'
```
