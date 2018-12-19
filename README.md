# Calanus
Albout-Boyer et al. 2016 Calanus GAMM

### To do

 + Start by focusing on stage 5 Calanus finmarchicus (Cfin_CV in the data spreadsheet)
 + Seasonal time series, ideally back to 1950
 + First stab: 
   1) build model using Caroline's data, but just SST and bathymetry, 
   2) use bathymetry and seasurface temperature (ERSST goes back to 1950) to project model
     - SST: https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.ersst.v5.html
     - BATHY: https://www.ngdc.noaa.gov/mgg/global/global.html
 + Nick will look up available models that give us more (e.g. bottom temperature, chlorophyll)
     - https://www.gfdl.noaa.gov/ocean-data-assimilation-model-output/
 + Add other Calanus data
   - *What are the units of the dataset we have?*
   - [ecomon](https://www.st.nmfs.noaa.gov/copepod/data/us-05101/)
   - [Jeff Runge](https://www.gmri.org/about-us/who-we-are/staff/jeffrey-runge-phd)

### Specifications

 + Bounding box
   -  left = -72, right = -42, bottom = 40, top = 70
 + Abundance units
   - Lehoux data... per m<sup>2</sup>
 
