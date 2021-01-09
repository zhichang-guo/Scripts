# Scripts
A repository for various scripts and tools
<!DOCTYPE html>
<html lang="en">
  <head></head>
  <body>
    <ol>
      <li>
      Find the nearest grid index with lat/lon stored in a netcdf file:
        <ul>
          <li>Command: python find_nearest_point.py -i input_file -l lon,lat</li>
          <li>Find the nearest point only over land points using "-g land"</li>
          <li>The default variable names for lon/lat/mask in the input_file are longitude/latitude/lsmask, otherwise specify them with "-v lon/lat/mask"</li>
          <li>An example: python find_nearest_point.py -i a.nc -l 60.10W,20.20S</li>
        </ul>
      </li>
    </ol>
  </body>
</html>
