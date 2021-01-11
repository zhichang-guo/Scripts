# Scripts
A repository for various scripts and tools
<!DOCTYPE html>
<html lang="en">
  <head></head>
  <body>
    <ol>
      Find the nearest grid index with given lat/lon from a pool of locations stored in a netcdf file:
        <li>1. For vector data, longitude/latitude are assumed to be functions of grid point.
                    This code will find the nearest grid point index with given lon/lat. 
                    python find_nearest_point.py -i file_name.nc -l 60W,20S</li>
             <li>2. For field data, lon/lat are assumed to be functions of i/j index respectively.
                    This code will find the nearest i/j index with given lon/lat.
                    python find_nearest_point.py -i file_name.nc -l 60W,20S -f field</li>
             <li>3. Find the nearest grid point only over land/ocean, land-sea mask is required. 
                    python find_nearest_point.py -i file_name.nc -l 60W,20S -g land</li>           
             <li>4. The default variable names are assumed to be longitude/latitude/lsmask. 
                    Otherwise, specify them with the flag -v.
                    python find_nearest_point.py -i file_name.nc -l 60W,20S -v lon/lat/mask</li>

          <li>Command: python find_nearest_point.py -i input_file -l lon,lat</li>
          <li>Find the nearest point only over land points using "-g land"</li>
          <li>The default variable names for lon/lat/mask in the input_file are longitude/latitude/lsmask, otherwise specify them with "-v lon/lat/mask"</li>
          <li>An example: python find_nearest_point.py -i a.nc -l 60.10W,20.20S</li>
       
      </li>
    </ol>
  </body>
</html>
