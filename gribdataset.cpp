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

#include <algorithm>

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

class GRIBAPIDataset : public GDALPamDataset
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
    FILE  *fp;
    grib_context *ctx;
    char  *pszProjection;
    double adfGeoTransform[6];
};

/************************************************************************/
/* ==================================================================== */
/*                            GRIBAPIRasterBand                             */
/* ==================================================================== */
/************************************************************************/

class GRIBAPIRasterBand : public GDALPamRasterBand
{
    friend class GRIBAPIDataset;

public:
    GRIBAPIRasterBand( GRIBAPIDataset *poDSIn, grib_handle *h );
    virtual ~GRIBAPIRasterBand();
    virtual CPLErr IReadBlock( int, int, void * );
    virtual const char *GetDescription() const;
    virtual double GetNoDataValue( int *pbSuccess = NULL );

private:
    grib_handle*  handle;

    int lastError;

    long nX, nY;

    double GetDouble ( const char *key, int *err ) const {
      double ret;
      *err = grib_get_double ( handle, key, &ret );
      return ret;
    }

    long GetLong ( const char *key, int *err ) const {
      long ret;
      *err = grib_get_long ( handle, key, &ret );
      return ret;
    }
};


/************************************************************************/
/*                           GRIBAPIRasterBand()                            */
/************************************************************************/

GRIBAPIRasterBand::GRIBAPIRasterBand( GRIBAPIDataset *poDSIn, grib_handle *handleIn ):
  handle(handleIn)
{
  lastError = GRIB_SUCCESS;
  poDS = poDSIn;

  nX = GetLong( "numberOfPointsAlongAParallel", &lastError );
  nY = GetLong( "numberOfPointsAlongAMeridian", &lastError );

  /*
    nBand = nBandIn;

    // Let user do -ot Float32 if needed for saving space, GRIBAPI contains
    // Float64 (though not fully utilized most of the time).
    eDataType = GDT_Float64;

    nBlockXSize = poDSIn->nRasterXSize;
    nBlockYSize = 1;

    const char* pszGribNormalizeUnits =
        CPLGetConfigOption("GRIB_NORMALIZE_UNITS", "YES");
    bool bMetricUnits = CPLTestBool(pszGribNormalizeUnits);

    SetMetadataItem( "GRIB_UNIT",
                     ConvertUnitInText(bMetricUnits, psInv->unitName) );
    SetMetadataItem( "GRIB_COMMENT",
                     ConvertUnitInText(bMetricUnits, psInv->comment) );
    SetMetadataItem( "GRIB_ELEMENT", psInv->element );
    SetMetadataItem( "GRIB_SHORT_NAME", psInv->shortFstLevel );
    SetMetadataItem( "GRIB_REF_TIME",
                     CPLString().Printf("%12.0f sec UTC", psInv->refTime ) );
    SetMetadataItem( "GRIB_VALID_TIME",
                     CPLString().Printf("%12.0f sec UTC", psInv->validTime ) );
    SetMetadataItem( "GRIB_FORECAST_SECONDS",
                     CPLString().Printf("%.0f sec", psInv->foreSec ) );
  */
}


/************************************************************************/
/*                         GetDescription()                             */
/************************************************************************/

const char * GRIBAPIRasterBand::GetDescription() const
{
    return GDALPamRasterBand::GetDescription();
}


/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr GRIBAPIRasterBand::IReadBlock( int /* nBlockXOff */,
                                   int nBlockYOff,
                                   void * pImage )

{
    return CE_None;
}

/************************************************************************/
/*                           GetNoDataValue()                           */
/************************************************************************/

double GRIBAPIRasterBand::GetNoDataValue( int *pbSuccess )
{
    *pbSuccess = FALSE;
    return 0;
}

/************************************************************************/
/*                           ~GRIBAPIRasterBand()                          */
/************************************************************************/

GRIBAPIRasterBand::~GRIBAPIRasterBand()
{
}

/************************************************************************/
/* ==================================================================== */
/*                              GRIBAPIDataset                             */
/* ==================================================================== */
/************************************************************************/

GRIBAPIDataset::GRIBAPIDataset():
  fp(0),
  ctx(0),
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
/*                            ~GRIBAPIDataset()                             */
/************************************************************************/

GRIBAPIDataset::~GRIBAPIDataset()

{
    FlushCache();
    CPLFree( pszProjection );
    if (ctx) {
      grib_context_delete(ctx);
    }
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

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *GRIBAPIDataset::Open( GDALOpenInfo * poOpenInfo )

{
    if( !Identify(poOpenInfo) )
        return NULL;

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

    poDS->fp = fopen(poOpenInfo->pszFilename, "r"); 

    if (!poDS->fp) {
        CPLDebug( "GRIB", "%s", strerror(errno) );
        CPLError( CE_Failure, CPLE_OpenFailed,
                  "Error (%d) opening file %s", errno, poOpenInfo->pszFilename);
        delete poDS;
        return NULL;
    }

    /*
    poDS->ctx = grib_context_new(grib_context_get_default());
    if (!poDS->ctx) {
        CPLError( CE_Failure, CPLE_OpenFailed,
                  "Error creating grib_context %s", poOpenInfo->pszFilename);
        delete poDS; // will close fp if not NULL
        return NULL;
    }
    */
    poDS->ctx = NULL; //FIXME

/* -------------------------------------------------------------------- */
/*      Create band objects.                                            */
/* -------------------------------------------------------------------- */
    int err;
    int bandNr = 1;
    GRIBAPIRasterBand *fstBand = NULL;

    grib_handle *h = grib_handle_new_from_file(poDS->ctx, poDS->fp, &err);
    if ( err != GRIB_SUCCESS ) goto open_error;
    fstBand = new GRIBAPIRasterBand ( poDS, h);
    if ( ( err = fstBand->lastError) != GRIB_SUCCESS ) {
      delete fstBand;
      goto open_error;
    }

    while ( ( h = grib_handle_new_from_file(poDS->ctx, poDS->fp, &err) ) != NULL ) {
      if ( err != GRIB_SUCCESS ) goto open_error;
      GRIBAPIRasterBand *band = new GRIBAPIRasterBand ( poDS, h );
      if ( ( err = band->lastError ) != GRIB_SUCCESS ) {
        if ( bandNr == 1 ) {
          delete fstBand;
        }
        goto open_error;
      }
      if ( bandNr == 1 ) {
        if ( fstBand->nX != band->nX || fstBand->nY != band->nY )  {
          CPLError( CE_Warning, CPLE_AppDefined,
                    "First message differs in size from the next, will skip it"
                    );
          delete fstBand;
          fstBand = band;
        } else {
          poDS->SetBand( bandNr++, fstBand );
        }
      }
      poDS->SetBand( bandNr++, band );
    }

/* -------------------------------------------------------------------- */
/*      Read common metadata from first valid message                         */
/* -------------------------------------------------------------------- */
    poDS->nRasterXSize = fstBand->nX;
    poDS->nRasterYSize = fstBand->nY;

/* -------------------------------------------------------------------- */
/*      Initialize any PAM information.                                 */
/* -------------------------------------------------------------------- */
    poDS->SetDescription( poOpenInfo->pszFilename );

    poDS->TryLoadXML();

/* -------------------------------------------------------------------- */
/*      Check for external overviews.                                   */
/* -------------------------------------------------------------------- */
    poDS->oOvManager.Initialize( poDS, poOpenInfo->pszFilename,
                                 poOpenInfo->GetSiblingFiles() );

    return poDS;

open_error:
  CPLError( CE_Failure, CPLE_OpenFailed,
            "Error (%s) opening file %s",
            grib_get_error_message(err),
            poOpenInfo->pszFilename 
          );
  delete poDS;
  return NULL;
}

/************************************************************************/
/*                            SetMetadata()                             */
/************************************************************************/

/*
void GRIBAPIDataset::SetGribMetaData(grib_MetaData* meta)
{
    nRasterXSize = meta->gds.Nx;
    nRasterYSize = meta->gds.Ny;

// --------------------------------------------------------------------
//      Image projection.                                              
// --------------------------------------------------------------------
    OGRSpatialReference oSRS;

    switch(meta->gds.projType)
    {
      case GS3_LATLON:
      case GS3_GAUSSIAN_LATLON:
          // No projection, only latlon system (geographic)
          break;
      case GS3_MERCATOR:
        oSRS.SetMercator(meta->gds.meshLat, meta->gds.orientLon,
                         1.0, 0.0, 0.0);
        break;
      case GS3_POLAR:
        oSRS.SetPS(meta->gds.meshLat, meta->gds.orientLon,
                   meta->gds.scaleLat1,
                   0.0, 0.0);
        break;
      case GS3_LAMBERT:
        oSRS.SetLCC(meta->gds.scaleLat1, meta->gds.scaleLat2,
                    meta->gds.meshLat, meta->gds.orientLon,
                    0.0, 0.0); // set projection
        break;

      case GS3_ORTHOGRAPHIC:

        // oSRS.SetOrthographic( 0.0, meta->gds.orientLon,
        //                       meta->gds.lon2, meta->gds.lat2);

        // oSRS.SetGEOS( meta->gds.orientLon, meta->gds.stretchFactor,
        //               meta->gds.lon2, meta->gds.lat2);

        // TODO: Hardcoded for now. How to parse the meta->gds section?
        oSRS.SetGEOS(  0, 35785831, 0, 0 );
        break;
      case GS3_EQUATOR_EQUIDIST:
        break;
      case GS3_AZIMUTH_RANGE:
        break;
    }

// --------------------------------------------------------------------
//      Earth model                                                   
// --------------------------------------------------------------------
    double a = meta->gds.majEarth * 1000.0; // in meters
    double b = meta->gds.minEarth * 1000.0;
    if( a == 0 && b == 0 )
    {
        a = 6377563.396;
        b = 6356256.910;
    }

    if (meta->gds.f_sphere)
    {
        oSRS.SetGeogCS( "Coordinate System imported from GRIB file",
                        NULL,
                        "Sphere",
                        a, 0.0 );
    }
    else
    {
        const double fInv = a / (a - b);
        oSRS.SetGeogCS( "Coordinate System imported from GRIB file",
                        NULL,
                        "Spheroid imported from GRIB file",
                        a, fInv );
    }

    OGRSpatialReference oLL; // construct the "geographic" part of oSRS
    oLL.CopyGeogCSFrom( &oSRS );

    double rMinX = 0.0;
    double rMaxY = 0.0;
    double rPixelSizeX = 0.0;
    double rPixelSizeY = 0.0;
    if (meta->gds.projType == GS3_ORTHOGRAPHIC)
    {
        // This is what should work, but it doesn't .. Dx seems to have an
        // inverse relation with pixel size.
        // rMinX = -meta->gds.Dx * (meta->gds.Nx / 2);
        // rMaxY = meta->gds.Dy * (meta->gds.Ny / 2);
        // Hardcoded for now, assumption: GEOS projection, full disc (like MSG).
        const double geosExtentInMeters = 11137496.552;
        rMinX = -(geosExtentInMeters / 2);
        rMaxY = geosExtentInMeters / 2;
        rPixelSizeX = geosExtentInMeters / meta->gds.Nx;
        rPixelSizeY = geosExtentInMeters / meta->gds.Ny;
    }
    else if( oSRS.IsProjected() )
    {
        // Longitude in degrees, to be transformed to meters (or degrees in
        // case of latlon).
        rMinX = meta->gds.lon1;
        // Latitude in degrees, to be transformed to meters.
        rMaxY = meta->gds.lat1;
        OGRCoordinateTransformation *poTransformLLtoSRS =
            OGRCreateCoordinateTransformation( &(oLL), &(oSRS) );
        // Transform it to meters.
        if( (poTransformLLtoSRS != NULL) &&
            poTransformLLtoSRS->Transform( 1, &rMinX, &rMaxY ))
        {
            if (meta->gds.scan == GRIBAPIBIT_2) // Y is minY, GDAL wants maxY
            {
                // -1 because we GDAL needs the coordinates of the centre of
                // the pixel.
                rMaxY += (meta->gds.Ny - 1) * meta->gds.Dy;
            }
            rPixelSizeX = meta->gds.Dx;
            rPixelSizeY = meta->gds.Dy;
        }
        else
        {
            rMinX = 0.0;
            rMaxY = 0.0;

            rPixelSizeX = 1.0;
            rPixelSizeY = -1.0;

            oSRS.Clear();

            CPLError( CE_Warning, CPLE_AppDefined,
                      "Unable to perform coordinate transformations, so the "
                      "correct projected geotransform could not be deduced "
                      "from the lat/long control points.  "
                      "Defaulting to ungeoreferenced." );
        }
        delete poTransformLLtoSRS;
    }
    else
    {
        // Longitude in degrees, to be transformed to meters (or degrees in
        // case of latlon).
        rMinX = meta->gds.lon1;
        // Latitude in degrees, to be transformed to meters.
        rMaxY = meta->gds.lat1;

        double rMinY = meta->gds.lat2;
        if (meta->gds.lat2 > rMaxY)
        {
          rMaxY = meta->gds.lat2;
          rMinY = meta->gds.lat1;
        }

        if( meta->gds.Nx == 1 )
          rPixelSizeX = meta->gds.Dx;
        else if (meta->gds.lon1 > meta->gds.lon2)
          rPixelSizeX =
              (360.0 - (meta->gds.lon1 - meta->gds.lon2)) / (meta->gds.Nx - 1);
        else
          rPixelSizeX = (meta->gds.lon2 - meta->gds.lon1) / (meta->gds.Nx - 1);

        if( meta->gds.Ny == 1 )
            rPixelSizeY = meta->gds.Dy;
        else
            rPixelSizeY = (rMaxY - rMinY) / (meta->gds.Ny - 1);

        // Do some sanity checks for cases that can't be handled by the above
        // pixel size corrections. GRIB1 has a minimum precision of 0.001
        // for latitudes and longitudes, so we'll allow a bit higher than that.
        if (rPixelSizeX < 0 || fabs(rPixelSizeX - meta->gds.Dx) > 0.002)
          rPixelSizeX = meta->gds.Dx;

        if (rPixelSizeY < 0 || fabs(rPixelSizeY - meta->gds.Dy) > 0.002)
          rPixelSizeY = meta->gds.Dy;
    }

    // http://gdal.org/gdal_datamodel.html :
    //   we need the top left corner of the top left pixel.
    //   At the moment we have the center of the pixel.
    rMinX-=rPixelSizeX/2;
    rMaxY+=rPixelSizeY/2;

    adfGeoTransform[0] = rMinX;
    adfGeoTransform[3] = rMaxY;
    adfGeoTransform[1] = rPixelSizeX;
    adfGeoTransform[5] = -rPixelSizeY;

    CPLFree( pszProjection );
    pszProjection = NULL;
    oSRS.exportToWkt( &(pszProjection) );
}
*/

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
    poDriver->SetMetadataItem( GDAL_DMD_EXTENSION, "grb" );
    //poDriver->SetMetadataItem( GDAL_DCAP_VIRTUALIO, "YES" );

    poDriver->pfnOpen = GRIBAPIDataset::Open;
    poDriver->pfnIdentify = GRIBAPIDataset::Identify;
    poDriver->pfnUnloadDriver = GDALDeregister_GRIBAPI;

    GetGDALDriverManager()->RegisterDriver( poDriver );
}
