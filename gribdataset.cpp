/******************************************************************************
 *
 * Project:  GRIB Driver
 * Purpose:  GDALDataset driver for GRIB using grib_api
 * Author:   Alberto Valverde, <alberto@meteogrid.com>
 *
 ******************************************************************************
 * Copyright (c) 2016, Alberto Valverde
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************
 *
 */
#include "grib_api.h"

#include "cpl_multiproc.h"
#include "gdal_frmts.h"
#include "gdal_pam.h"
#include "ogr_spatialref.h"

#include <assert.h>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#define MAX_KEY_LEN  255
#define MAX_VAL_LEN  1024

CPL_CVSID("$Id$");

CPL_C_START
void GDALRegister_GRIBAPI();
CPL_C_END

/************************************************************************/
/* ==================================================================== */
/*                              GRIBAPIDataset                             */
/* ==================================================================== */
/************************************************************************/

class GRIBAPIRasterBand;

typedef std::vector<GRIBAPIRasterBand*> band_vector;

class GRIBAPIDataset : public GDALDataset
{
    friend class GRIBAPIRasterBand;

  public:
                GRIBAPIDataset();
                ~GRIBAPIDataset();

    static GDALDataset *Open( GDALOpenInfo * );
    static int          Identify( GDALOpenInfo * );

    CPLErr      GetGeoTransform( double * padfTransform );
    const char *GetProjectionRef();

  private:
    CPLErr LoadMetaData(const band_vector &bands);
    CPLErr LoadSubdatasets(const band_vector &bands);

    GRIBAPIRasterBand *GetGribBand ( const int i );

    FILE  *fp;
    char  *pszProjection;
    double adfGeoTransform[6];
};

/************************************************************************/
/* ==================================================================== */
/*                            GRIBAPIRasterBand                             */
/* ==================================================================== */
/************************************************************************/

class GRIBAPIRasterBand : public GDALRasterBand
{
    friend class GRIBAPIDataset;

public:
  GRIBAPIRasterBand( GRIBAPIDataset *, int, grib_handle * );
  virtual ~GRIBAPIRasterBand();
  virtual CPLErr IReadBlock( int, int, void * );
  virtual const char *GetDescription() const;
  virtual double GetNoDataValue( int *pbSuccess = NULL );


private:
  grib_handle*  handle;

  long nX, nY;

  CPLErr GetGeoTransform (double (&padfTransform)[6]) const;

  CPLErr LoadMetaData( const char* name_space );
  CPLErr LoadMetaData( ) { return LoadMetaData ( NULL ); }

  double GetDouble ( const char *key, int *err ) const {
    double ret;
    *err = grib_get_double ( handle, key, &ret );
    if ( *err != GRIB_SUCCESS ) {
      CPLError( CE_Failure, CPLE_AppDefined,
                "GetDouble(%s) failed: %s\n", key,
                grib_get_error_message(*err) );
    }
    return ret;
  }

  long GetLong ( const char *key, int *err ) const {
    long ret;
    *err = grib_get_long ( handle, key, &ret );
    if ( *err != GRIB_SUCCESS ) {
      CPLError( CE_Failure, CPLE_AppDefined,
                "GetLong(%s) failed: %s\n", key,
                grib_get_error_message(*err) );
    }
    return ret;
  }

  long GetLong ( const char *key ) const {
    int err;
    return GetLong ( key, &err);
  }

  int GetString ( const char *key, char *value, size_t *vlen ) const {
    int err;
    err = grib_get_string ( handle, key, value, vlen );
    if ( err != GRIB_SUCCESS ) {
      CPLError( CE_Failure, CPLE_AppDefined,
                "GetString(%s) failed: %s\n", key,
                grib_get_error_message(err) );
    }
    return err;
  }

  const std::string GetString ( const char *key ) const {
    size_t len = MAX_VAL_LEN;
    char buffer[MAX_VAL_LEN];
    if ( GetString ( key, buffer, &len ) != GRIB_SUCCESS ) {
      return std::string("<INVALID_KEY>");
    } else {
      return std::string(buffer);
    }
  }

  std::string GetIdentifier () const {
    if ( GetLong( "level" ) ) {
      return GetString( "level" ) + GetString( "indicatorOfParameter" );
    } else {
      return GetString( "indicatorOfParameter" );
    }
  }

  static bool allSameSize (const band_vector& bands)
  {
    band_vector::const_iterator it = bands.begin();
    GRIBAPIRasterBand *fst = *(it++);
    for ( ; it < bands.end(); it++ ) {
      if ( (*it)->nX != fst->nX || (*it)->nY != fst->nY ) {
        return false;
      }
    }
    return true;
  }

};


/************************************************************************/
/*                           GRIBAPIRasterBand()                            */
/************************************************************************/

GRIBAPIRasterBand::GRIBAPIRasterBand(
    GRIBAPIDataset *poDSIn, int nBandIn, grib_handle *handleIn ):
  handle(handleIn)
{
  poDS = poDSIn;
  nBand = nBandIn;
  nX = nBlockXSize = GetLong( "Ni" );
  nY = nBlockYSize = GetLong( "Nj" );
  eDataType = GDT_Float64;
}

/************************************************************************/
/*                         LoadMetaData()                               */
/************************************************************************/
CPLErr GRIBAPIRasterBand::LoadMetaData( const char *name_space )
{
  unsigned long filter_flags =
    CPLTestBool( CPLGetConfigOption( "GRIBAPI_ALL_KEYS", "OFF" ) ) ?
    GRIB_KEYS_ITERATOR_ALL_KEYS : GRIB_KEYS_ITERATOR_SKIP_COMPUTED;

  grib_keys_iterator* kiter =
    grib_keys_iterator_new( handle, filter_flags, name_space );
  if ( !kiter ) {
    CPLError( CE_Failure, CPLE_OpenFailed,
              "GRIBAPI: grib_keys_iterator_new returned NULL\n" );
    return CE_Failure;
  }
  while( grib_keys_iterator_next( kiter ) ) {
    char value[MAX_VAL_LEN];
    const char* name = grib_keys_iterator_get_name(kiter);
    if (!strncmp(name, "value", strlen("value")) ||
        !strncmp(name, "pv", strlen("pv"))) continue;

    size_t vlen = MAX_VAL_LEN;
    if ( GRIB_SUCCESS == GetString ( name, value, &vlen ) ) {
      std::ostringstream key;
      key << "GRIB_";
      if ( name_space ) {
        std::string ns( name_space );
        std::transform(ns.begin(), ns.end(), ns.begin(), ::toupper);
        key << ns << "_";
      }
      key << name;
      SetMetadataItem ( key.str().c_str(), value );
    }
  }
  grib_keys_iterator_delete ( kiter );
  return CE_None;
}


/************************************************************************/
/*                         GetDescription()                             */
/************************************************************************/

const char * GRIBAPIRasterBand::GetDescription() const
{
    return GDALRasterBand::GetDescription();
}

static inline double degToRad(double a) {
  return a * M_PI / 180;
}

/************************************************************************/
/*                         GetGeoTransform()                             */
/************************************************************************/
CPLErr GRIBAPIRasterBand::GetGeoTransform (double (&padfTransform)[6]) const
{
  int err;
  const double degRot = GetDouble( "angleOfRotationInDegrees", &err );
  if (err) return CE_Failure;
  if ( degRot != 0 ) {
    CPLError( CE_Warning, CPLE_NotSupported,
              "The GRIBAPI driver does not support yet "
              "datasets with rotation.\n" );
  }
  //
  // Longitude in degrees, to be transformed to meters (or degrees in
  // case of latlon).
  const double lon1 = GetDouble( "longitudeOfFirstGridPointInDegrees", &err );
  if (err) return CE_Failure;

  const double lon2 = GetDouble( "longitudeOfLastGridPointInDegrees", &err );
  if (err) return CE_Failure;

  const double lat1 = GetDouble( "latitudeOfFirstGridPointInDegrees", &err );
  if (err) return CE_Failure;

  const double lat2 = GetDouble( "latitudeOfLastGridPointInDegrees", &err );
  if (err) return CE_Failure;

  const double dlon = GetDouble( "iDirectionIncrementInDegrees", &err );
  if (err) return CE_Failure;

  const double dlat = GetDouble( "jDirectionIncrementInDegrees", &err );
  if (err) return CE_Failure;

  const bool iScansPos = GetLong( "iScansPositively", &err );
  if (err) return CE_Failure;

  const bool jScansPos = GetLong( "jScansPositively", &err );
  if (err) return CE_Failure;

  double rPixelSizeX = 0.0;

  if ( nX == 1 )
    rPixelSizeX = dlon;
  else if (lon1 > lon2)
    rPixelSizeX = (360.0 - (lon1 - lon2)) / (nX - 1);
  else
    rPixelSizeX = (lon2 - lon1) / (nX - 1);

  double rPixelSizeY = 0.0;

  if( nY == 1 )
      rPixelSizeY = dlat;
  else if (lat1 > lat2)
      rPixelSizeY = (lat1 - lat2) / (nY - 1);
  else
      rPixelSizeY = (lat2 - lat1) / (nY - 1);

  // Do some sanity checks for cases that can't be handled by the above
  // pixel size corrections. GRIB1 has a minimum precision of 0.001
  // for latitudes and longitudes, so we'll allow a bit higher than that.
  if (rPixelSizeX < 0 || fabs(rPixelSizeX - dlon) > 0.002)
    rPixelSizeX = dlon;

  if (!iScansPos)
    rPixelSizeX *= -1;

  if (rPixelSizeY < 0 || fabs(rPixelSizeY - dlat) > 0.002)
    rPixelSizeY = dlat;

  if (!jScansPos)
    rPixelSizeY *= -1;

  // move to center of pixel

  padfTransform[0] = lon1 - rPixelSizeX/2;
  padfTransform[1] = rPixelSizeX;
  padfTransform[2] = 0;
  padfTransform[3] = lat1 - rPixelSizeY/2;
  padfTransform[4] = 0;
  padfTransform[5] = rPixelSizeY;

  return CE_None;
}


/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr GRIBAPIRasterBand::IReadBlock( int /* nBlockXOff */,
                                      int /* nBlockYOff */,
                                      void * pImage )
{
  size_t size;

  int err;
  if ( ( err = grib_get_size ( handle, "values", &size ) ) != GRIB_SUCCESS )  {
    CPLError( CE_Failure, CPLE_AppDefined,
              "GRIBAPI: Could not get value array size: %s.\n",
              grib_get_error_message(err));
    return CE_Failure;
  }

  if ( size != nX*nY ) {
    CPLError( CE_Failure, CPLE_AppDefined,
              "GRIBAPI: unexpected value array size.\n" );
    return CE_Failure;
  }

  double *buff = static_cast<double *>(pImage);
  if ( ( err = grib_get_double_array ( handle, "values", buff, &size ) ) != GRIB_SUCCESS )  {
    CPLError( CE_Failure, CPLE_AppDefined,
              "GRIBAPI: Could get valye array: %s.\n",
              grib_get_error_message(err));
    return CE_Failure;
  }
    return CE_None;
}

/************************************************************************/
/*                           GetNoDataValue()                           */
/************************************************************************/

double GRIBAPIRasterBand::GetNoDataValue( int *pbSuccess )
{
    int err;
    double result = GetDouble("missingValue", &err);
    if ( pbSuccess )
      *pbSuccess = err == GRIB_SUCCESS ? TRUE : FALSE;
    return result;
}

/************************************************************************/
/*                           ~GRIBAPIRasterBand()                          */
/************************************************************************/

GRIBAPIRasterBand::~GRIBAPIRasterBand()
{
  grib_handle_delete (handle);
}

/************************************************************************/
/*                              GRIBAPIDataset                             */
/************************************************************************/

GRIBAPIDataset::GRIBAPIDataset():
  fp(0),
  pszProjection(0)
{
  adfGeoTransform[0] = 0.0;
  adfGeoTransform[1] = 1.0;
  adfGeoTransform[2] = 0.0;
  adfGeoTransform[3] = 0.0;
  adfGeoTransform[4] = 0.0;
  adfGeoTransform[5] = 1.0;
}

/************************************************************************/
/*                              GetGribBand                             */
/************************************************************************/
GRIBAPIRasterBand *GRIBAPIDataset::GetGribBand ( const int i )
{
  return static_cast<GRIBAPIRasterBand*>( GetRasterBand(i) );
}

/************************************************************************/
/*                              LoadMetaData                            */
/************************************************************************/
CPLErr GRIBAPIDataset::LoadMetaData(const band_vector& bands)
{
/* -------------------------------------------------------------------- */
/*      Read common metadata from first valid message                         */
/* -------------------------------------------------------------------- */
  CPLErr err = CE_None;
  if (bands.size()) {

    for (band_vector::const_iterator it=bands.begin(); it<bands.end(); it++) {
      if ( ( err = (*it)->LoadMetaData() ) != CE_None ) {
        return err;
      }
    }
    nRasterXSize = bands[0]->nX;
    nRasterYSize = bands[0]->nY;
    if ( ( err = bands[0]->GetGeoTransform( adfGeoTransform ) ) != CE_None )
      return err;
  }

  CPLFree( pszProjection );
  pszProjection = NULL;
  OGRSpatialReference oSRS;
  oSRS.importFromEPSG(4326);
  oSRS.exportToWkt( &(pszProjection) );

  return err;
}

CPLErr GRIBAPIDataset::LoadSubdatasets(const band_vector& bands)
{
  band_vector::const_iterator it;
  int i;
  for (i=0, it=bands.begin(); it<bands.end(); i++, it++) {
    GRIBAPIRasterBand *band = *it;
    std::ostringstream name_key;
    name_key << "SUBDATASET_" << i + 1 << "_NAME" ;
    std::ostringstream name_value;
    name_value << "GRIBAPI:"
               << band->GetIdentifier()
               << ":"
               << GetDescription();
    SetMetadataItem ( name_key.str().c_str(),
                      name_value.str().c_str(),
                      "SUBDATASETS" );

    std::ostringstream desc_key;
    desc_key << "SUBDATASET_" << i + 1 << "_DESC" ;
    std::ostringstream desc_value;
    desc_value
      << band->GetString("parameterName")
      << " (" << band->GetString("dataDate")
              << band->GetString("dataTime") << " +"
              << band->GetString("startStep")
              << band->GetString("stepUnits")
      << ") of "
      << GetDescription();
    SetMetadataItem ( desc_key.str().c_str(),
                      desc_value.str().c_str(),
                      "SUBDATASETS" );
  }
  return CE_None;
}

/************************************************************************/
/*                            ~GRIBAPIDataset()                             */
/************************************************************************/

GRIBAPIDataset::~GRIBAPIDataset()

{
    FlushCache();
    CPLFree( pszProjection );
    if (fp) {
      fclose(fp);
    }
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr GRIBAPIDataset::GetGeoTransform( double * padfTransform )

{
    memcpy( padfTransform,  adfGeoTransform, sizeof(double) * 6 );
    return CE_None;
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char *GRIBAPIDataset::GetProjectionRef()

{
    return pszProjection;
}

/************************************************************************/
/*                            Identify()                                */
/************************************************************************/

int GRIBAPIDataset::Identify( GDALOpenInfo * poOpenInfo )
{
  const char *pszFilename = poOpenInfo->pszFilename;
  if( STARTS_WITH_CI(pszFilename, "GRIBAPI:") ) {
    int i;
    for ( i = 8; i<strlen(pszFilename); i++) {
      if (pszFilename[i] == ':') {
        GDALOpenInfo oOpenInfo( pszFilename + i + 1, poOpenInfo->eAccess );
        return Identify(&oOpenInfo);
      }
    }
    return FALSE;
  }

  if (poOpenInfo->nHeaderBytes < 8)
    return FALSE;

/* -------------------------------------------------------------------- */
/*      Does a part of what ReadSECT0() but in a thread-safe way.       */
/* -------------------------------------------------------------------- */
    for(int i=0;i<poOpenInfo->nHeaderBytes-3;i++)
    {
        if (STARTS_WITH_CI((const char*)poOpenInfo->pabyHeader + i, "GRIB") ||
            STARTS_WITH_CI((const char*)poOpenInfo->pabyHeader + i, "TDLP"))
            return TRUE;
    }

    return FALSE;
}

static inline void deleteBands(band_vector& bands)
{
  band_vector::const_iterator it;
  for ( it = bands.begin(); it < bands.end(); it++ ) delete *it;
  bands.erase(bands.begin(), bands.end());
}


/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *GRIBAPIDataset::Open( GDALOpenInfo * poOpenInfo )

{
  if( !Identify(poOpenInfo) )
    return NULL;

  const char *pszFilename = poOpenInfo->pszFilename;
  char *subdataset = NULL;
  long subdatasetLong = 0;
  if( STARTS_WITH_CI(pszFilename, "GRIBAPI:") ) {
    int i;
    for ( i = 8; i<strlen(pszFilename); i++) {
      if (pszFilename[i] == ':') {
        subdataset = strndup(pszFilename+8, i-8);
        subdatasetLong = atoi(subdataset);
        pszFilename += i + 1;
        break;
      }
    }
  }

/* -------------------------------------------------------------------- */
/*      Confirm the requested access is supported.                      */
/* -------------------------------------------------------------------- */
    if( poOpenInfo->eAccess == GA_Update )
    {
        CPLError( CE_Failure, CPLE_NotSupported,
                  "The GRIBAPI driver does not support update access to existing"
                  " datasets.\n" );
        return NULL;
    }
/* -------------------------------------------------------------------- */
/*      Create a corresponding GDALDataset.                             */
/* -------------------------------------------------------------------- */
    GRIBAPIDataset *poDS = new GRIBAPIDataset();

    poDS->SetDescription( pszFilename );

    poDS->fp = fopen(pszFilename, "r"); 

    if (!poDS->fp) {
        CPLDebug( "GRIB", "%s", strerror(errno) );
        CPLError( CE_Failure, CPLE_OpenFailed,
                  "Error (%d) opening file %s\n",
                  errno, pszFilename );
        delete poDS;
        return NULL;
    }

/* -------------------------------------------------------------------- */
/*      Create band objects.                                            */
/* -------------------------------------------------------------------- */
    int err;
    grib_handle *h = grib_handle_new_from_file(NULL, poDS->fp, &err);

    if ( err != GRIB_SUCCESS ) {
      CPLError( CE_Failure, CPLE_OpenFailed,
                "Error (%s) opening file %s\n",
                grib_get_error_message(err),
                poOpenInfo->pszFilename 
              );
      delete poDS;
      return NULL;
    }

    band_vector bands;
    // Create bands from all messages if no subdataset requested or only
    // for the requested one
    while ( ( h = grib_handle_new_from_file(NULL, poDS->fp, &err) ) != NULL ) {
      if ( err != GRIB_SUCCESS ) {
        CPLError( CE_Failure, CPLE_OpenFailed,
                  "Error (%s) opening file %s\n",
                  grib_get_error_message(err),
                  poOpenInfo->pszFilename 
                );
        delete poDS;
        return NULL;
      }
      GRIBAPIRasterBand *band = new GRIBAPIRasterBand ( poDS, 0, h );
      if (subdataset) {
        if ( !strcmp( subdataset, band->GetIdentifier().c_str() ) ) {
          bands.push_back(band);
        } else {
          delete band;
        }
      } else {
          bands.push_back(band);
      }
    }

    // free it now that we no longer use it. Make sure not to use the pointer
    // from now on!
    if (subdataset) {
      free(subdataset);
    }
    
    if ( bands.size() == 0 ) {
      CPLError( CE_Warning, CPLE_AppDefined,
                "GRIBAPI: Could not find any GRIB message\n"
                );
      delete poDS;
      return NULL;
      
    }

    // INM sends a dummy first message, skip it
    if ( !subdataset && !GRIBAPIRasterBand::allSameSize ( bands ) ) {
      CPLError( CE_Warning, CPLE_AppDefined,
                "GRIBAPI: First message differs in size from the next. "
                "I will skip it.\n"
                );
      delete bands[0];
      bands.erase( bands.begin() );
      if ( !GRIBAPIRasterBand::allSameSize ( bands ) ) {
        CPLError( CE_Warning, CPLE_AppDefined,
                  "GRIBAPI: Second message differs too. I can't handle this\n"
            );
        deleteBands( bands );
        delete poDS;
        return NULL;
      }
    }

    if ( CE_None != poDS->LoadMetaData(bands) ) {
      CPLError( CE_Failure, CPLE_OpenFailed,
                "Could not load metadata from %s\n",
                poOpenInfo->pszFilename );
      deleteBands( bands );
      delete poDS;
      return NULL;
    }

    if ( !subdataset && CE_None != poDS->LoadSubdatasets(bands) ) {
      CPLError( CE_Failure, CPLE_OpenFailed,
                "Could not query subdatasets from  %s\n",
                poOpenInfo->pszFilename );
      deleteBands( bands );
      delete poDS;
      return NULL;
    }

    if (subdataset) {
      // Add the bands from the vector to the dataset only if a subdataset
      // was requested
      int i;
      band_vector::const_iterator it;
      for ( it = bands.begin(), i = 1 ; it < bands.end() ; it++, i++) {
        (*it)->nBand = i;
        poDS->SetBand( i, *it );
      }
    }

    return poDS;
}


/************************************************************************/
/*                       GDALDeregister_GRIBAPI()                          */
/************************************************************************/

static void GDALDeregister_GRIBAPI(GDALDriver* )
{
}

/************************************************************************/
/*                         GDALRegister_GRIBAPI()                          */
/************************************************************************/

void GDALRegister_GRIBAPI()

{
    if( GDALGetDriverByName( "GRIBAPI" ) != NULL )
        return;

    GDALDriver *poDriver = new GDALDriver();

    poDriver->SetDescription( "GRIBAPI" );
    poDriver->SetMetadataItem( GDAL_DCAP_RASTER, "YES" );
    poDriver->SetMetadataItem( GDAL_DMD_LONGNAME, "GRIdded Binary (.grb)" );
    poDriver->SetMetadataItem( GDAL_DMD_HELPTOPIC, "frmt_grib2.html" );
    poDriver->SetMetadataItem( GDAL_DMD_SUBDATASETS, "YES" );
    poDriver->SetMetadataItem( GDAL_DMD_EXTENSION, "grb" );
    //poDriver->SetMetadataItem( GDAL_DCAP_VIRTUALIO, "YES" );

    poDriver->pfnOpen = GRIBAPIDataset::Open;
    poDriver->pfnIdentify = GRIBAPIDataset::Identify;
    poDriver->pfnUnloadDriver = GDALDeregister_GRIBAPI;

    GetGDALDriverManager()->RegisterDriver( poDriver );
}
