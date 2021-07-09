import xarray as xr
import numpy as np


class InfluenceFunction(object):
    def __init__(self, filename, delimiter='_', verbose=False,
                 latname='lat', lonname='lon', zname=None):
        # Should add multi-file options
        self.VERBOSE = verbose
        self.latname = latname
        self.lonname = lonname
        self.zname = zname
        self.dataset = xr.open_dataset(filename)
        self.delimiter = delimiter
        self.collect_meta()

    def collect_meta(self):
        self.sources = []
        self.species = []
        self.pathways = []

        for var in self.dataset.variables:
            try:
                source, species, pathway = self.varname_to_meta(var)
                if source not in self.sources:
                    self.sources.append(source)
                if species not in self.species:
                    self.species.append(species)
                if pathway not in self.pathways:
                    self.pathways.append(pathway)
            except ValueError:
                if self.VERBOSE:
                    print(f'Variable not elligible: {var}')
                else:
                    pass

        if self.VERBOSE:
            print(f'Source Types: {self.sources}')
            print(f'Chemical Species: {self.species}')
            print(f'Exposure Pathways: {self.pathways}')

    def get_latitudes(self):
        return self.dataset[self.latname].values

    def get_longitudes(self):
        return self.dataset[self.lonname].values

    def get_lonslats(self):
        return self.get_longitudes(), self.get_latitudes()

    def get_altitudes(self):
        if self.zname is None:
            print("No altitude axis")
            return None
        else:
            return self.dataset[self.zname].values

    def varname_to_meta(self, varname):
        source, species, pathway = varname.split(self.delimiter)
        return source, species, pathway

    def meta_to_varname(self, source, species, pathway):
        return f'{source}{self.delimiter}{species}{self.delimiter}{pathway}'

    def get_IF(self, source, species, pathway):
        varname = self.meta_to_varname(source, species, pathway)
        return np.squeeze(self.dataset[varname].values)

    def get_species_exposure(self, source_dict, species, pathway):
        exposure = self.get_IF(
            self.sources[0], self.species[0], self.pathways[0]) * 0.0
        for source, s in source_dict.items():
            exposure += self.get_IF(source, species, pathway) * s
        return exposure
