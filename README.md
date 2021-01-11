# Scripts
A repository for various scripts and tools
<!DOCTYPE html>
<html lang="en">
  <head></head>
  <body>
    <ol>
      <h3>Find the nearest grid index from locations stored in a netcdf file with provided lat/lon:</h3>
        <ul><li>For vector data, longitude/latitude are assumed to be functions of grid point.
                This code will find the nearest grid point index with given lon/lat.<br> 
                python find_nearest_point.py -i file_name.nc -l 60W,20S</li>
            <li>For field data, lon/lat are assumed to be functions of i/j index respectively.
                This code will find the nearest i/j index with given lon/lat.<br>
                python find_nearest_point.py -i file_name.nc -l 60W,20S -f field</li>
            <li>Find the nearest grid point only over land/ocean, land-sea mask is required.<br> 
                python find_nearest_point.py -i file_name.nc -l 60W,20S -g land</li>           
            <li>The default variable names are assumed to be longitude/latitude/lsmask. 
                Otherwise, specify them with the flag -v.<br>
                python find_nearest_point.py -i file_name.nc -l 60W,20S -v lon/lat/mask</li>
        </ul>
      <h3>Plot time series of variables stored in netcdf files, examples:</h3>
        <ul><li>python plot_timeseries_nc.py -i file_name.nc -v vname1,vname2 <br>
                1D variable is assumed to be time-dependent. The code draws time series of 2 variables.
                2D(vector,time) variable, it draws time series of 2 variables averaged over all points.
                3D(lon,lat,time) variable. The code draws domain average time series of 2 variables.</li>
            <li>python plot_timeseries_nc.py -i data_name.nc -v vname1-vname2 -t 2,-2 -factor 86400<br>
                Draw time series of difference, factored by 86400, between vname1 and vname2 for the
                time segment from the second record to the second last record.</li>
            <li>python plot_timeseries_nc.py -i data_name.nc -g geo_name -v vname -lon lon -lat lat<br>
                The code will find the grid point which is closest to the given location(lon/lat) and 
                draw the time series. The variable vname is assumed to be stored in data_name.nc while
                the geographic information is stored in geo_name. If the files are located in the same
                directory, the directory path just needs to be specified in one file name.</li>
            <li>python plot_timeseries_nc.py -i data_name.nc -v vname -lon lon1,lon2 -lat lat1,lat2<br>
                The code will average the variable over the domain (lon1,lon2),(lat1,lat2) and draw
                the time series.</li>
            <li>python plot_timeseries_nc.py -i data_name.nc -v vname -p point_index -l dotted -c red<br>
                Draw time series of the variable at the location with point_index. The data are drawn
                with the dotted red line.</li>
            <li>python plot_timeseries_nc.py -i data_name.nc -v vname -p point_index -s no<br>
                The code draws the time series and output it to a png file.</li>
        </ul>
    </ol>
  </body>
</html>
