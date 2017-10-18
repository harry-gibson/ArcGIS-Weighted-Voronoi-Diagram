cimport cython, openmp
from libc.math cimport sqrt, sin, cos, atan, M_PI, tan, atan2, fabs, asin, ceil
from cython.parallel cimport parallel, prange
import numpy as np

cdef:
    int MAX_ITERATIONS = 200
    double CONVERGENCE_THRESHOLD = 1e-12
    double TO_RADIANS = M_PI / 180.0
    double V_a = 6378137
    double V_f = 1 / 298.257223563
    double V_b = 6356752.314245
    double AVG_EARTH_RADIUS = 6371

cdef class weighted_voronoi:

    cdef:
        Py_ssize_t rasterWidth, rasterHeight, nPoints
        double[::1] lats, lons, pointWts, pointLats, pointLons
        int[::1] pointIds

    def __cinit__(self, double[::1] lats, double[::1] lons,
                  int[::1] pointIds, double[::1] pointWts,
                  double[::1] pointLats, double[::1] pointLons):
        self.rasterHeight = len(lats)
        self.rasterWidth = len(lons)
        self.nPoints = len(pointIds)
        assert((len(pointWts)==self.nPoints) and (len(pointLats) == self.nPoints) and
               (len(pointLons)==self.nPoints))
        self.lats = lats
        self.lons = lons
        self.pointIds = pointIds
        self.pointWts = pointWts
        self.pointLats = pointLats
        self.pointLons = pointLons


    cdef float haversine(self, float y1, float x1, float y2, float x2) nogil:
        """ Calculate the great-circle distance bewteen two points on the Earth surface.
    
        :input: four floats, containing the latitude and longitude of each point
        in decimal degrees (lat1, lon1, lat2, lon2)
    
        Example: haversine(45.7597, 4.8422, 48.8567, 2.3508)
    
        :output: Returns the distance bewteen the two points.
        The default unit is kilometers. Miles can be returned
        if the ``miles`` parameter is set to True.
    
        """

        cdef:
            double xr1, xr2, yr1, yr2, d, h, xDist, yDist

        xr1 = x1 * TO_RADIANS
        yr1 = y1 * TO_RADIANS
        xr2 = x2 * TO_RADIANS
        yr2 = y2 * TO_RADIANS


        # calculate haversine
        yDist = yr2 - yr1
        xDist = xr2 - xr1
        d = sin(yDist * 0.5) ** 2 + cos(yr1) * cos(yr2) * sin(xDist * 0.5) ** 2
        h = 2 * AVG_EARTH_RADIUS * asin(sqrt(d))
        return h  # in kilometers

    @cython.cdivision(True)
    cdef float vincenty(self, float y1, float x1, float y2, float x2) nogil:
        cdef:
            Py_ssize_t iteration
            double U1, U2, L, Lambda, sinU1, cosU1, sinU2, cosU2
            double sinLambda, cosLambda, sinSigma, cosSigma, sigma
            double sinAlpha, cosSqAlpha, cos2SigmaM, LambdaPrev
            double C
            double uSq, A, B, deltaSigma, s

        if x1 == x2 and y1 == y2:
            return 0.0

        U1 = atan((1 - V_f) * tan(y1 * TO_RADIANS))
        U2 = atan((1 - V_f) * tan(y2 * TO_RADIANS))
        L = (x2 - x1) * TO_RADIANS
        Lambda = L

        sinU1 = sin(U1)
        cosU1 = cos(U1)
        sinU2 = sin(U2)
        cosU2 = cos(U2)

        for iteration in range (MAX_ITERATIONS):
            sinLambda = sin(Lambda)
            cosLambda = cos(Lambda)
            sinSigma = sqrt((cosU2 * sinLambda) ** 2 +
                                 (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
            if sinSigma == 0:
                return 0.0  # coincident points
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
            sigma = atan2(sinSigma, cosSigma)
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
            cosSqAlpha = 1 - sinAlpha ** 2
            if cosSqAlpha != 0:
                cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
            else:
                cos2SigmaM = 0
            C = V_f / 16 * cosSqAlpha * (4 + V_f * (4 - 3 * cosSqAlpha))
            LambdaPrev = Lambda
            Lambda = L + (1 - C) * V_f * sinAlpha * (sigma + C * sinSigma *
                                                   (cos2SigmaM + C * cosSigma *
                                                    (-1 + 2 * cos2SigmaM ** 2)))
            if fabs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
                break  # successful convergence
        else:
            return -1
        uSq = cosSqAlpha * (V_a ** 2 - V_b ** 2) / (V_b**2)
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
        deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                     (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                     (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
        s = V_b * A * (sigma - deltaSigma)

        s /= 1000  # meters to kilometers

        return s


    @cython.boundscheck(False)
    @cython.cdivision(True)
    cpdef int[:,::1] calculate(self, float power = 1, unsigned char useEllipsoid=1):
        cdef:
            Py_ssize_t rasterWidth, rasterHeight, x, y, nPoints, pt
            int[:,::1] assignments
            int assignment = -1
            float minScore, curScore, inf
            int fivePercentRows
            double[::1] lats, lons, pointWts, pointLats, pointLons
            int[::1] pointIds

        rasterWidth = self.rasterWidth
        rasterHeight = self.rasterHeight
        lats = self.lats
        lons = self.lons
        pointIds = self.pointIds
        pointWts = self.pointWts
        pointLats = self.pointLats
        pointLons = self.pointLons

        nPoints = self.nPoints
        #fivePercentRows = int(ceil(rasterHeight / 20) + 1)

        assignments = np.empty(shape=(rasterHeight, rasterWidth),
                                    dtype=np.int32)

        assignments = np.empty(shape=(rasterHeight,rasterWidth), dtype=np.int32)
        inf = np.inf
        with nogil, cython.wraparound(False), parallel():
            for y in range(rasterHeight):
                for x in prange(rasterWidth):
                    assignment = -1
                    minScore = inf
                    curScore = inf
                    for pt in range(nPoints):
                        # todo allow alternate distance weighting e.g.
                        # f(d) = 100 + d^(3/2); the scale is 400 by 600.
                        if useEllipsoid:
                            curScore =  self.vincenty(pointLats[pt], pointLons[pt],
                                                            lats[y], lons[x]) / pointWts[pt]
                        else:
                            curScore = self.haversine(pointLats[pt], pointLons[pt],
                                                            lats[y], lons[x]) / pointWts[pt]

                        if curScore < minScore:
                            minScore = curScore
                            assignment = pointIds[pt]
                    assignments[y, x] = assignment
            return assignments

