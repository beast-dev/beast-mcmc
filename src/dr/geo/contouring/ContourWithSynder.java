package dr.geo.contouring;

import dr.geo.KernelDensityEstimator2D;

import java.util.*;
import java.awt.geom.Point2D;


/**
 * @author Marc A. Suchard
 */
public class ContourWithSynder extends KernelDensityEstimator2D  {

    public ContourWithSynder(double[] x, double[] y, double[] h, int n, double[] lims) {
        super(x, y, h, n, lims);
    }

    public ContourWithSynder(double[] x, double[] y) {
        super(x, y);
    }

    public ContourWithSynder(double[] x, double[] y, int n) {
        super(x, y, n);
    }

    public ContourPath[] getContourPaths(double hpdValue) {

        if (contourPaths == null) {

            double thresholdDensity = findLevelCorrespondingToMass(hpdValue);

            SnyderContour contourPlot = new SnyderContour(getXGrid().length,getYGrid().length);
            contourPlot.setDeltas(getXGrid()[1]-getXGrid()[0],getYGrid()[1]-getYGrid()[0] );
            contourPlot.setOffsets(getXGrid()[0],getYGrid()[0]);

            List<LinkedList<Point2D>> allPaths = new ArrayList<LinkedList<Point2D>>();
            contourPlot.ContourKernel(getKDE(),allPaths,thresholdDensity);

            contourPaths = new ContourPath[allPaths.size()];
            for(int i=0; i<allPaths.size(); i++) {
                LinkedList<Point2D> path = allPaths.get(i);
                int len = path.size();
                double[] x = new double[len];
                double[] y = new double[len];
                for(int j=0; j<len; j++) {
                    Point2D pt = path.get(j);
                    x[j] = pt.getX();
                    y[j] = pt.getY();
                }
                contourPaths[i] = new ContourPath(new ContourAttrib(thresholdDensity),1,x,y);
            }
        }

        return contourPaths;
    }

    private ContourPath[] contourPaths = null;

}
