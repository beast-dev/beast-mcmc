package dr.geo;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.stats.DiscreteStatistics;

import java.util.Arrays;

/**
 * KernelDensityEstimator2D creates a bi-variate kernel density smoother for data
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class KernelDensityEstimator2D {

//    kde2d =
//    function (x, y, h, n = 25, lims = c(range(x), range(y)))
//    {
//        nx <- length(x)
//        if (length(y) != nx)
//            stop("data vectors must be the same length")
//        if (any(!is.finite(x)) || any(!is.finite(y)))
//            stop("missing or infinite values in the data are not allowed")
//        if (any(!is.finite(lims)))
//            stop("only finite values are allowed in 'lims'")
//        gx <- seq.int(lims[1], lims[2], length.out = n)
//        gy <- seq.int(lims[3], lims[4], length.out = n)
//        if (missing(h))
//            h <- c(bandwidth.nrd(x), bandwidth.nrd(y))
//        h <- h/4
//        ax <- outer(gx, x, "-")/h[1]
//        ay <- outer(gy, y, "-")/h[2]
//        z <- matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n,
//            nx))/(nx * h[1] * h[2])
//        return(list(x = gx, y = gy, z = z))
//    }

    /*
     * @param x x-coordinates of observations
     * @param y y-coordinates of observations
     * @param h bi-variate smoothing bandwidths
     * @param n smoothed grid size
     * @param lims bi-variate min/max for grid
     */
    public KernelDensityEstimator2D(double[] x, double[] y, double[] h, int n, double[] lims) {
        this.x = x;
        this.y = y;
        if (x.length != y.length)
            throw new RuntimeException("data vectors must be the same length");

        this.nx = x.length;

        if (n <= 0)
            throw new RuntimeException("must have a positive number of grid points");
        this.n = n;

        if (lims != null)
            this.lims = lims;
        else
            setupLims();

        if (h != null)
            this.h = h;
        else
            setupH();

        doKDE2D();
    }

    public KernelDensityEstimator2D(double[] x, double[] y) {
        this(x,y,null,25,null);
    }

    public KernelDensityEstimator2D(double[] x, double[] y, int n) {
        this(x,y,new double[]{1.0,1.0},n,null);
    }

    public void doKDE2D() {
        gx = makeSequence(lims[0], lims[1], n);
        gy = makeSequence(lims[2], lims[3], n);
        double[][] ax = outerMinusScaled(gx, x, h[0]);
        double[][] ay = outerMinusScaled(gy, y, h[1]);
        normalize(ax);
        normalize(ay);
        z = new double[n][n];
        double scale = nx * h[0] * h[1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double value = 0;
                for (int k = 0; k < nx; k++) {
                    value += ax[i][k] * ay[j][k];
                }
                z[i][j] = value / scale;
            }
        }
    }

    public double findLevelCorrespondingToMass(double probabilityMass) {
        double level = 0;
        double[] sz = new double[n*n];
        double[] c1 = new double[n*n];
        for(int i=0; i<n; i++)
            System.arraycopy(z[i],0,sz,i*n,n);
        Arrays.sort(sz);
        final double dx = gx[1] - gx[0];
        final double dy = gy[1] - gy[0];
        final double dxdy = dx * dy;
        c1[0] = sz[0] * dxdy;
        final double criticalValue = 1.0 - probabilityMass;
        if (criticalValue < c1[0] || criticalValue >= 1.0)
                throw new RuntimeException();
        // do linearInterpolation on density (y) as function of cummulative sum (x)
        for(int i=1; i<n*n; i++) {
            c1[i] = sz[i] * dxdy + c1[i-1];
            if (c1[i] > criticalValue) { // first largest point
                final double diffC1 = c1[i] - c1[i-1];
                final double diffSz = sz[i] - sz[i-1];
                level = sz[i] - (c1[i]-criticalValue) / diffC1 * diffSz;
                break;
            }
        }
        return level;
    }

    public double[][] getKDE() {
        return z;
    }

    public double[] getXGrid() {
        return gx;
    }

    public double[] getYGrid() {
        return gy;
    }

    public void normalize(double[][] X) {
        for (int i = 0; i < X.length; i++) {
            for (int j = 0; j < X[0].length; j++)
                X[i][j] = NormalDistribution.pdf(X[i][j], 0, 1);
        }
    }

    public double[][] outerMinusScaled(double[] X, double[] Y, double scale) {
        double[][] A = new double[X.length][Y.length];
        for (int indexX = 0; indexX < X.length; indexX++) {
            for (int indexY = 0; indexY < Y.length; indexY++)
                A[indexX][indexY] = (X[indexX] - Y[indexY]) / scale;
        }
        return A;
    }

    public double[] makeSequence(double start, double end, int length) {
        double[] seq = new double[length];
        double by = (end - start) / (length - 1);
        double value = start;
        for (int i = 0; i < length; i++, value += by) {
            seq[i] = value;
        }
        return seq;
    }

    private void setupLims() {
        lims = new double[4];
        lims[0] = DiscreteStatistics.min(x);
        lims[1] = DiscreteStatistics.max(x);
        lims[2] = DiscreteStatistics.min(y);
        lims[3] = DiscreteStatistics.max(y);
    }

    private void setupH() {
        h = new double[2];
        h[0] = bandwidthNRD(x) / 4;
        h[1] = bandwidthNRD(y) / 4;
    }


//   bandwidth.nrd =
//   function (x)
//   {
//       r <- quantile(x, c(0.25, 0.75))
//       h <- (r[2] - r[1])/1.34
//       4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)

    //   }
    public double bandwidthNRD(double[] in) {

        DoubleArrayList inList = new DoubleArrayList(in.length);
        for (double d : in)
            inList.add(d);
        inList.sort();

        final double h = (Descriptive.quantile(inList, 0.75) - Descriptive.quantile(inList, 0.25)) / 1.34;

        return 4 * 1.06 *
                Math.min(Math.sqrt(DiscreteStatistics.variance(in)), h) *
                Math.pow(in.length, -0.2);
    }

//    public double[][] constructContour(double height)
//    {
//
//        if (height <= minHeight || height >= maxHeight) {
//            return false;
//        }
//
//        double z1, z2, z3, z4;
//        double x1, x2, y1, y2;
//        double[] topX = new double[1];
//        double[] rightY = new double[1];
//        double[] bottomX = new double[1];
//        double[] leftY = new double[1];
//        double[] middleX = new double[1];
//        double[] middleY = new double[1];
//        boolean top, right, bottom, left, middle;
//
//        int count;
//        for (int x = 1; x < columnCount; x++) {
//            x1 = mXScale->getVal(x-1);
//            x2 = mXScale->getVal(x);
//
//            for (int y = 1; y < rowCount; y++) {
//                y1=mYScale->getVal(y-1);
//                y2=mYScale->getVal(y);
//
//                z1=mSurface->getElement(y-1, x-1);
//                z2=mSurface->getElement(y-1, x);
//                z3=mSurface->getElement(y, x);
//                z4=mSurface->getElement(y, x-1);
//
//                count = 0;
//                if (top = intersectSide(height, x1, z4, x2, z3, topX)) {
//                    count++;
//                }
//
//                if (left = intersectSide(height, y1, z1, y2, z4, leftY)) {
//                    count++;
//                }
//
//                if (middle = intersectSide(height, y1, z1, y2, z3, middleY)) {
//                    intersectSide(height, x1, z1, x2, z3, middleX);
//                    count++;
//                }
//
//                if (count==1) {
//                    throw new RuntimeException("error constructing contour");
//                }
//
//                if (top) {
//                    tmpXData->add(topX);
//                    tmpYData->add(y2);
//                }
//                if (left) {
//                    tmpXData->add(x1);
//                    tmpYData->add(leftY);
//                }
//                if (middle) {
//                    tmpXData->add(middleX);
//                    tmpYData->add(middleY);
//                }
//
//                count = 0;
//                if (bottom = intersectSide(height, x1, z1, x2, z2, bottomX)) {
//                    count++;
//                }
//
//                if (right = intersectSide(height, y1, z2, y2, z3, rightY)) {
//                    count++;
//                }
//
//                if (middle) {
//                    count++;
//                }
//
//                if (count==1) {
//                    throw new RuntimeException("error constructing contour");
//                }
//
//                if (bottom) {
//                    tmpXData->add(bottomX);
//                    tmpYData->add(y1);
//                }
//                if (right) {
//                    tmpXData->add(x2);
//                    tmpYData->add(rightY);
//                }
//                if (middle) {
//                    tmpXData->add(middleX);
//                    tmpYData->add(middleY);
//                }
//            }
//        }
//
//        connectContour(tmpXData, tmpYData, outXData, outYData);
//
//        return true;
//    }
//
//    void connectContour(MEDataColumn *inXData, MEDataColumn *inYData, MEDataColumn *outXData, MEDataColumn *outYData)
//    {
//        // Find a segment that is not connected to another
//        boolean foundEnd, foundMatch;
//        int i, j, n=inXData->getN();
//
//        i=0;
//        foundEnd=false;
//        while (!foundEnd && i<n) {
//            j=i+2;
//            foundMatch=false;
//            while (!foundMatch && j<n) {
//                if (inXData->getVal(i)==inXData->getVal(j) &&
//                        inYData->getVal(i)==inYData->getVal(j))
//                    foundMatch=true;
//
//                j+=2;
//            }
//            if (!foundMatch)
//                foundEnd=true;
//            else
//                i+=2;
//        }
//
//        // If foundEnd is false then the contour is closed - start anywhere (i.e. at 0)
//        if (!foundEnd) {
//            i=0;
//            do {
//                outXData->add(inXData->getVal(i));
//                outYData->add(inYData->getVal(i));
//
//                j=0;
//                foundMatch=false;
//                while (!foundMatch && j<n) {
//                    if (j!=i) {
//                        if (inXData->getVal(i)==inXData->getVal(j) &&
//                                inYData->getVal(i)==inYData->getVal(j))
//                            foundMatch=true;
//                    }
//                    j+=2;
//                }
//
//                ThrowIf(!foundMatch);
//
//                i=j;
//            } while (i!=0);
//
//            outXData->add(inXData->getVal(0));
//            outYData->add(inYData->getVal(0));
//
//        } else {
//            foundEnd=false;
//            do {
//                outXData->add(inXData->getVal(i));
//                outYData->add(inYData->getVal(i));
//
//                j=0;
//                foundMatch=false;
//                while (!foundMatch && j<n) {
//                    if (j!=i) {
//                        if (inXData->getVal(i)==inXData->getVal(j) &&
//                                inYData->getVal(i)==inYData->getVal(j))
//                            foundMatch=true;
//                    }
//                    j+=2;
//                }
//
//                if (foundMatch) {
//                    i=j;
//                } else {
//                    outXData->add(inXData->getVal(i+1));
//                    outYData->add(inYData->getVal(i+1));
//
//                    foundEnd=true;
//                }
//            } while (!foundEnd);
//        }
//    }

    private boolean intersectSide(double z0, double u1, double z1, double u2, double z2, double []v)
    {
        double du, dz, ddz;

        if (z1<z2 && z0>z1 && z0<z2) {
            du = u2 - u1;
            dz = z2 - z1;
            ddz = z0 - z1;
            v[0] = u1 + ((du * ddz) / dz);
//		v=u1 + (((u2 - u1) * (z0 - z1)) / (z2 - z1));
            return true;
        } else if (z2<z1 && z0>z2 && z0<z1) {
            du = u2 - u1;
            dz = z1 - z2;
            ddz = z0 - z2;
            v[0] = u2 - ((du * ddz) / dz);
//		v=u1 + (((u2 - u1) * (z0 - z2)) / (z1 - z2));
            return true;
        }

        return false;
    }


    public static void main(String[] arg) {

        double[] x = {3.4, 1.2, 5.6, 2.2, 3.1};
        double[] y = {1.0, 2.0, 1.0, 2.0, 1.0};

        KernelDensityEstimator2D kde = new KernelDensityEstimator2D(x, y, 4);

        System.out.println(new Vector(kde.getXGrid()));
        System.out.println(new Vector(kde.getYGrid()));
        System.out.println(new Matrix(kde.getKDE()));
        System.exit(-1);

    }

    private double[] x; // x coordinates
    private double[] y; // y coordinates
    private double[] h; // h[0] x-bandwidth, h[1] y-bandwidth
    private int n; // grid size
    private double[] lims; // x,y limits
    private int nx; // length of vectors
    private double[] gx; // x-grid points
    private double[] gy; // y-grid points
    private double[][] z; // KDE estimate;

}
