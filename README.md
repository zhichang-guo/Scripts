# Scripts
A repository for various scripts and tools
<!DOCTYPE html>
<html lang="en">
  <head></head>
  <body>
    <ol>
      <li>
      Find the nearest grid index with given lat/lon from a pool of locations stored in a netcdf file:
        <ul><li>For vector data, longitude/latitude are assumed to be functions of grid point.
                This code will find the nearest grid point index with given lon/lat. 
                python find_nearest_point.py -i file_name.nc -l 60W,20S</li>
            <li>For field data, lon/lat are assumed to be functions of i/j index respectively.
                This code will find the nearest i/j index with given lon/lat.
                python find_nearest_point.py -i file_name.nc -l 60W,20S -f field</li>
            <li>Find the nearest grid point only over land/ocean, land-sea mask is required. 
                python find_nearest_point.py -i file_name.nc -l 60W,20S -g land</li>           
            <li>The default variable names are assumed to be longitude/latitude/lsmask. 
                Otherwise, specify them with the flag -v.
                python find_nearest_point.py -i file_name.nc -l 60W,20S -v lon/lat/mask</li>
        </ul>
      </li>
    </ol>
  </body>
</html>
