from time import time
# needs to be before the arcpy import as that screws up the paths
from weighted_voronoi import weighted_voronoi
import arcpy
import numpy as np

startTime = time()

pointSource = arcpy.GetParameterAsText(0)
weightField = arcpy.GetParameterAsText(1)
vectorOutput = arcpy.GetParameterAsText(2)
bufferDegrees = float(arcpy.GetParameterAsText(3))
cellSize = float(arcpy.GetParameterAsText(4))
distMethod = arcpy.GetParameterAsText(5)
smoothPolys = True if arcpy.GetParameterAsText(6) == "true" else False

# Split out the destination file into path and name
vectorOutput = vectorOutput.replace("\\", "/")
vectorOutputPath = vectorOutput[:vectorOutput.rindex("/")]
vectorOutputName = vectorOutput[vectorOutput.rindex("/")+1:]

# Get the OID field name
oidField = arcpy.Describe(pointSource).OIDFieldName


def cleanup_scratch(refs):
    """
    Cleanup intermediate files that may be in scratch
    """
    # Cleanup intermediate files that may be around
    for r in refs:
        arcpy.Delete_management(r)


def cell_center_lat(row, initlat, ysize, numrows):
    """
    The latitude for the cell centers in a row
    """
    return initlat + (numrows - row - 1) * ysize


def cell_center_lng(col, initlng, xsize, numcols):
    """
    The longitude for the cell centers in a column
    """
    return initlng + (col) * xsize


# Get the spatial reference of the source (and dest)
spatialRef = arcpy.Describe(pointSource).spatialReference

# Split out the destination file into path and name
vectorOutputPath = vectorOutput[:vectorOutput.rindex("/")]
vectorOutputName = vectorOutput[vectorOutput.rindex("/")+1:]

# Intermediate features/rasters
vectorIntermediate = arcpy.env.scratchGDB + "/DVWPoly_vect_int"
rasterIntermediate = arcpy.env.scratchGDB + "/DVWPoly_rast_int"
rasterVals = arcpy.env.scratchGDB + "/DVWPoly_rast_vals"
rasterIds = arcpy.env.scratchGDB + "/DVWPoly_rast_ids"
extentBuffered = arcpy.env.scratchGDB + "/DVWPoly_extent_buffered"

# References to anything that might be in scratch
scratchRefs = [vectorIntermediate, rasterIntermediate, rasterVals,
               rasterIds, extentBuffered]

# Cleanup anything that may be leftover in scratch
cleanup_scratch(scratchRefs)

# Ensure that a scratch GDB is set
if arcpy.env.scratchGDB is None:
    arcpy.AddError("No scratch geodatabase is set")
    cleanup_scratch(scratchRefs)
    raise arcpy.ExecuteError

# Ensure that the point source is in a GCS
if spatialRef.type != "Geographic":
    arcpy.AddError("Must use a geographic coordinate system")
    cleanup_scratch(scratchRefs)
    raise arcpy.ExecuteError

# Get the extent of the point source and create as a polygon
sourceExtent = arcpy.Describe(pointSource).extent
extentPoints = arcpy.Array()
extentPoints.add(sourceExtent.lowerLeft)
extentPoints.add(sourceExtent.upperLeft)
extentPoints.add(sourceExtent.upperRight)
extentPoints.add(sourceExtent.lowerRight)
extentPoints.add(sourceExtent.lowerLeft)
sourceExtentPoly = arcpy.Polygon(extentPoints)

# Buffer the point source extent by bufferDegrees and update
# processing extent environment variable
arcpy.Buffer_analysis(sourceExtentPoly, extentBuffered,
                      bufferDegrees, "OUTSIDE_ONLY")
initExtent = arcpy.env.extent
arcpy.env.extent = arcpy.Describe(extentBuffered).extent

# Create a raster of weight values from weightField
arcpy.FeatureToRaster_conversion(
    in_features=pointSource,
    field=weightField,
    out_raster=rasterVals,
    cell_size=cellSize)
valRaster = arcpy.Raster(rasterVals)

# Create a raster of OIDs
arcpy.FeatureToRaster_conversion(
    in_features=pointSource,
    field=oidField,
    out_raster=rasterIds,
    cell_size=cellSize)
oidRaster = arcpy.Raster(rasterIds)

# Return to the initial processing extent
arcpy.env.extent = initExtent

# Create value and OID numpy arrays from rasters
arcpy.AddMessage("Creating value and OID numpy arrays")
valArray = arcpy.RasterToNumPyArray(valRaster, nodata_to_value=0)
oidArray = arcpy.RasterToNumPyArray(oidRaster, nodata_to_value=-1)


lowerLeft = arcpy.Point(oidRaster.extent.XMin, oidRaster.extent.YMin)
llCellX = oidRaster.extent.XMin + 0.5 * oidRaster.meanCellWidth
llCellY = oidRaster.extent.YMin + 0.5 * oidRaster.meanCellHeight
cellSizeY = oidRaster.meanCellHeight
cellSizeX = oidRaster.meanCellWidth
rasterHeight = int(oidRaster.height)
rasterWidth = int(oidRaster.width)

arcpy.AddMessage("Rasters height: %d" % rasterHeight)
arcpy.AddMessage("Rasters width: %d" % rasterWidth)

del valRaster
del oidRaster

# Arrays of the cell center latitudes for each row and longitudes for
# each column
arcpy.AddMessage("Determine lat for each row, lng for each column")
rowLats = np.array(
                    [cell_center_lat(x, llCellY, cellSizeY, rasterHeight)
                     for x
                     in range(rasterHeight)])

colLngs = np.array(
                    [cell_center_lng(y, llCellX, cellSizeX, rasterWidth)
                     for y
                     in range(rasterWidth)])

# Get the row/col refs for point locations in the raster and build a
# list of points including OID and coordinates
arcpy.AddMessage("Building up points to check")
pointLocations = zip(*map(list, np.where(oidArray != -1)))
#checkPoints = [
#    {
#        "cell": l,
#        "coords": (rowLats[l[0]], colLngs[l[1]]),
#        "oid": oidArray[l],
#        "val": valArray[l]
#    }
#    for l in pointLocations
#]
oidToUse = oidArray[oidArray!=-1]
valToUse = valArray[oidArray!=-1]
latToUse = np.asarray([rowLats[l[0]] for l in pointLocations])
lonToUse = np.asarray([colLngs[l[1]] for l in pointLocations])

del oidArray
del valArray

# the cython has been compiled for int32 specifically
assignmentsDType = np.int32

doEllipsoid = 0 if distMethod == "Haversine" else 1

# do thinking
arcpy.AddMessage("Beginning distance allocation")
vCalc = weighted_voronoi(rowLats, colLngs, oidToUse, valToUse, latToUse, lonToUse)
assignments = np.asarray(vCalc.calculate(useEllipsoid=doEllipsoid))

# Save the array back to a raster
arcpy.AddMessage("Assignments completed, saving to raster")
tempRast = arcpy.NumPyArrayToRaster(assignments, lowerLeft,
                                   cellSizeX, cellSizeY)
tempRast.save(rasterIntermediate)
del tempRast

# Define the raster projection
arcpy.DefineProjection_management(rasterIntermediate, spatialRef)

# Create an intermediate polygon feature class in scratch because
# AlterField_management doesn't appear to work with the output shapefile.
arcpy.AddMessage("Polygonizing raster")
arcpy.RasterToPolygon_conversion(
    in_raster=rasterIntermediate,
    out_polygon_features=vectorIntermediate,
    simplify=("SIMPLIFY" if smoothPolys else "NO_SIMPLIFY"),
    raster_field="Value")

# Delete the kinda-duplicative Id field
arcpy.AddMessage("Cleaning up polygon fields")
arcpy.DeleteField_management(vectorIntermediate, "Id")

# Rename GRIDCODE field to align with the original OID fieldname, e.g. ORIG_FID
arcpy.AlterField_management(vectorIntermediate,
                            "GRIDCODE",
                            str("ORIG_" + oidField))

# Copy the intermediate vector output to its final destination
arcpy.AddMessage("Saving final polygon output to {}".format(vectorOutput))
arcpy.FeatureClassToFeatureClass_conversion(
    in_features=vectorIntermediate,
    out_path=vectorOutputPath,
    out_name=vectorOutputName
  )

# Cleanup intermediate files
cleanup_scratch(scratchRefs)

