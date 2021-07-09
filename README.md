# source-receptor-tools

[(See included PFAS example)](pfas_example.ipynb)

## Requirements
 - xarray
 - numpy

## Metadata Structure

Based on tree structure:

`Source Type -> Chemical Species -> Exposure Pathway -> Spatiotemporal Information`

Relies on NetCDF/xarray for handling the spatiotemporal dimensions.

Implemented with \*.nc files using the Variable convention:
```python
name = f'{SOURCETYPE}_{CHEMICALSPECIES}_{EXPOSUREPATHWAY}'
```
