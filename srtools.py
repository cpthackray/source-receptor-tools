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
        """Get latitude values of associated dataset.

        Returns:
        latitudes [array]
        """
        return self.dataset[self.latname].values

    def get_longitudes(self):
        """Get longitude values of associated dataset.

        Returns:
        longitudes [array]
        """
        return self.dataset[self.lonname].values

    def get_lonslats(self):
        """Get longitudes and latitudes of associated dataset.

        Returns:
        ( longitudes [array], latitudes [array] )
        """
        return self.get_longitudes(), self.get_latitudes()

    def get_altitudes(self):
        """Get altitudes of associated dataset (if present).

        Returns:
        altitudes [array or None]
        """
        if self.zname is None:
            print("No altitude axis")
            return None
        else:
            return self.dataset[self.zname].values

    def varname_to_meta(self, varname):
        """Convert file's variable name to required metadata set.

        Arguments:
        varname [string]

        Returns:
        (source, species, pathway) [strings]
        """
        source, species, pathway = varname.split(self.delimiter)
        return source, species, pathway

    def meta_to_varname(self, source, species, pathway):
        """Convert metadata set to file's variable name.

        Arguments:
        source [string]
        species [string]
        pathway [string]

        Returns:
        [string]
        """
        return f'{source}{self.delimiter}{species}{self.delimiter}{pathway}'

    def get_IF(self, source, species, pathway):
        """Retrive influence functions for given metadata set.

        Arguments:
        source [string]
        species [string]
        pathway [string]

        Returns:
        [Ntimes x Nlats x Nlons x Nalts array]
        """
        varname = self.meta_to_varname(source, species, pathway)
        return np.squeeze(self.dataset[varname].values)

    def get_species_exposure(self, source_dict, species, pathway):
        """Exposure for species and pathway associated with defined source magnitudes.

        Arguments:
        source_dict [dict] - contains sourcename:float pairs defining source magnitudes
        species [string]
        pathway [string]

        Returns:
        [Ntimes x Nlats x Nlons x Natls array]"""
        exposure = self.get_IF(
            self.sources[0], self.species[0], self.pathways[0]) * 0.0
        for source, s in source_dict.items():
            exposure += self.get_IF(source, species, pathway) * s
        return exposure
