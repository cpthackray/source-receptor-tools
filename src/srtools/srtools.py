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

    def lonlat_to_ij(self, lon, lat):
        lons, lats = self.get_lonslats()
        j = np.argmin(abs(lon - lons))
        i = np.argmin(abs(lat - lats))
        return i, j

    def exposure_at_locations(self, source_dict, species, pathway, locations):
        """Get exposure for specific lon/lat pairs."""
        exposure_map = self.get_species_exposure(source_dict=source_dict, species=species, pathway=pathway)
        exposures = np.zeros(len(locations))
        for iloc, (lon, lat) in enumerate(locations):
            i, j = self.lonlat_to_ij(lon, lat)
            exposures[iloc] = exposure_map[i,j]
        return exposures

    def fill_dataframe(self, df, source_dict, species_pathway_list,
                        latname = 'latitude', lonname = 'longitude'):
        """Fill dataframe that has locations with listed (species, pathway) pairs."""
        locations = [(lon, lat) for lon, lat in zip(df[lonname].values, df[latname].values)]
        for spc, pathway in species_pathway_list:
            exposures = self.exposure_at_locations(source_dict=source_dict, species=spc, 
                                                    pathway=pathway, locations=locations)
            df[f'{spc}_{pathway}'] = exposures
        return df

    def fill_netcdf(self, source_dict, species_pathway_list, filename):
        dv = {}
        for spc, pw in species_pathway_list:
            arr = self.get_species_exposure(source_dict=source_dict, species=spc, pathway=pw)
            dv[f'{spc}_{pw}'] = (('time','lat','lon'),np.expand_dims(arr,0))
        ds = xr.Dataset(data_vars=dv,
                    coords={'lat':(('lat'),self.get_latitudes()),
                            'lon':(('lon'),self.get_longitudes()),
                            'time':(('time'),[0.0])})
        ds.to_netcdf(filename)
