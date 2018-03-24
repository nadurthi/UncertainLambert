function Y=convert_ecef2longlatalt_coordinates(Y,constants)


[long,lat,Alt]=cspice_recgeo(Y',constants.radii(1),constants.f);
