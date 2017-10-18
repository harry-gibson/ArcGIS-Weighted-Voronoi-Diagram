try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import  Extension

from Cython.Build import cythonize
import os

# https://github.com/cython/cython/wiki/PackageHierarchy
def scandir(dir, files=[]):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            # get as a dot-separated path excluding the extensio
            files.append((path.replace(os.path.sep, ".")[:-4]))
        elif os.path.isdir(path):
            # recursive scan
            scandir(path, files)
    return files


def makeExtension(extName):
    extPath = extName.replace(".", os.path.sep)+".pyx"
    return Extension(
        extName,
        [extPath],
        include_dirs=["."],
        extra_compile_args=['/openmp', '-O3'],
        extra_link_args=['/openmp']
    )

#voronoiExtNames = scandir("weighted_voronoi_c")
extensions = [Extension(
          "weighted_voronoi",
          ["weighted_voronoi.pyx"],
          extra_compile_args=['/openmp'],
          extra_link_args=['/openmp'],
          include_dirs=['.']
   )]

#allextNames = spatialextNames + temporalextNames
#extensions = [makeExtension(name) for name in voronoiExtNames]

setup(
    name = "Weighted Voronoi allocation",
    description = "Raster-based weighted voronoi allocation",
    author = "Harry Gibson",
    author_email = "harry.s.gibson@gmail.com",
    long_description='''
    Cythonised implementation of the core algorithm from:
    https://github.com/wiringa/ArcGIS-Weighted-Voronoi-Diagram
    This is for weighted voronoi polygon generation, implemented using a raster-based 
    algorithm that calculates at each cell location the lowest-cost input point, where 
    a lower cost is a function of both (decreased) distance and (increased) value of a 
    given weight field. 
    ''',
    #packages=["weighted_voronoi_c", "weighted_voronoi_c.geodesic_dist"],
    #ext_modules = cythonize(extensions),
    ext_modules=extensions,
    #py_modules=['aggregation_values',
    #            'spatial.spatial_aggregation_runner',
    #            'temporal.temporal_aggregation_runner']
)