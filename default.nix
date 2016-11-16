{ stdenv, gdal, grib-api, cmake }:

stdenv.mkDerivation rec {
  version = "1.0";
  shortname = "gdal_gribapi";
  name = "${shortname}-${version}";

  src = ./.;

  nativeBuildInputs = [ cmake ];
  cmakeFlags = [];

  buildInputs = [ gdal grib-api ];

  doCheck = true;

  checkPhase = ''
    export GDAL_SKIP=GRIB
    export GDAL_DRIVER_PATH="$(pwd)"
    gdalinfo --formats | grep GRIBAPI
    '';

  meta = {
    description = "GDAL plugin to read GRIB files with grib_api";
    homepage = https://github.com/meteogrid/gdal_gribapi;
    license = stdenv.lib.licenses.bsd3;
  };
}
