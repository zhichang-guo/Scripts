# Scripts
A repository for scripts: <a href="#nearest">Find the nearest point</a>, <a href="#timeseries">Time series plot</a>, <a href="#snapshot">Geographical plot</a>
<!DOCTYPE html>
<html lang="en">
  <head></head>
  <body>
    <ol>
      <h3><a id="user-content-nearest" href="#nearest">Given longitude/latitude, find the nearest grid point from locations stored in a netcdf file</a>:</h3>
        <ul><li>For vector data, longitude/latitude are assumed to be functions of grid point.
                This code will find the nearest grid point index with lon/lat provided.<br> 
                <pre><code>python find_nearest_point.py -i file_name.nc -l 60W,20S</code></pre></li>
            <li>For field data, lon/lat are assumed to be functions of i/j index respectively.
                This code will find the nearest i/j index with given lon/lat.<br>
                <pre><code>python find_nearest_point.py -i file_name.nc -l 60W,20S -f field</code></pre></li>
            <li>Find the nearest grid point only over land or ocean. It is noted that land-sea mask is required.<br> 
                <pre><code>python find_nearest_point.py -i file_name.nc -l 60W,20S -g land</code></pre></li>           
            <li>The default variable names are longitude/latitude/lsmask. Otherwise, specify them with the flag -v.<br>
                <pre><code>python find_nearest_point.py -i file_name.nc -l 60W,20S -v lon,lat,mask</code></pre></li>
        </ul>
      <h3><a id="user-content-timeseries" href="#timeseries">Plot time series of variables stored in netcdf files, examples</a>:</h3>
        <ul><li><pre><code>python plot_timeseries_nc.py -i file_name.nc -v vname1,vname2 </code></pre>
                For 1D (time) variable, the code draws time series of 2 variables;
                For 2D(vector,time) variable, it draws time series of 2 variables averaged over all points;
                For 3D(lon,lat,time) variable, it draws domain average time series of 2 variables.</li>
            <li><pre><code>python plot_timeseries_nc.py -i data_name.nc -v vname1-vname2 -t 2,-3 -factor 86400</code></pre>
                It draws time series of difference, multiplied by 86400, between vname1 and vname2 for the
                time segment from the second record to the third last record.</li>
            <li><pre><code>python plot_timeseries_nc.py -i data.nc -g geo.nc -v vname -lon lon -lat lat</code></pre>
                The code will find the grid point which is closest to the given location(lon/lat) and 
                draw the time series. The variable vname is assumed to be stored in data.nc while
                the geographic information is stored in geo.nc. If the files are located in the same
                directory, the directory path just needs to be specified in one file.</li>
            <li><pre><code>python plot_timeseries_nc.py -i data_name.nc -v vname -lon lon1,lon2 -lat lat1,lat2</code></pre>
                The code will average the variable over the domain (lon1,lon2),(lat1,lat2) and draw
                the time series.</li>
            <li><pre><code>python plot_timeseries_nc.py -i data_name.nc -v vname -p point_index -l dotted -c red</code></pre>
                Draw time series of the variable at the location with point_index. The data are drawn
                with the dotted red line.</li>
            <li><pre><code>python plot_timeseries_nc.py -i data_name.nc -v vname -p point_index -s no</code></pre>
                The code draws the time series and outputs it to a png file.</li>
        </ul>
      <h3><a id="user-content-snapshot" href="#snapshot">Create a plot on a map with irregularly spaced data stored in netcdf files, examples</a>:</h3>
        <ul><li><pre><code>python plot_snapshot_nc.py -i file_name.nc -v vname</code></pre> 
                1D variable is assumed to be stationary vector data. Also assume lon/lat and the variable
                are stored in the same file. For 2D(vector,time) and 3D(lon,lat,time) variable, the 
                code creates a geographical plot with data averaged over the time domain.</li>
            <li><pre><code>python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -f field</code></pre> 
                The variable vname is assumed to be stored in data.nc while the geographic information is stored in
                geo.nc. If the files are located in the same directory, the directory path just needs to be
                specified in one file. With the flag "-f field", 2D variable is assumed to be stationary field data (lon,lat)
                instead of time-varying vectors (vector,time).</li>
            <li><pre><code>python plot_snapshot_nc.py -i data_name.nc -v vname1-vname2 -t 2,-3 -factor 86400</code></pre>
                Make a graphical plot of difference between vname1 and vname2. The data are multiplied by 86400 and
                averaged over the time segment from the second record to the third last record.</li>
            <li><pre><code>python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -lon lon1,lon2 -lat lat1,lat2</code></pre>
                The code will make a geographical plot over a specific spatial domain.</li>
            <li><pre><code>python plot_snapshot_nc.py -i data.nc -g geo.nc -v vname -s no -t 128</code></pre>
                The code makes a geographical plot for the time step 128 and output it to a png file.</li>
        </ul>
    </ol>
  </body>
</html>
