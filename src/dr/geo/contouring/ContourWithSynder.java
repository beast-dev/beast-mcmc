/*
 * ContourWithSynder.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.geo.contouring;

import dr.geo.KernelDensityEstimator2D;

import java.util.*;
import java.awt.geom.Point2D;


/**
 * @author Marc A. Suchard
 */
public class ContourWithSynder extends KernelDensityEstimator2D  {

    public ContourWithSynder(final double[] x, final double[] y, final double[] h, final int n, final double[] lims) {
        super(x, y, h, n, lims);
    }

    public ContourWithSynder(final double[] x, final double[] y, boolean bandwidthLimit) {
        super(x, y, bandwidthLimit);
    }

    public ContourWithSynder(final double[] x, final double[] y) {
        super(x, y);
    }

    public ContourWithSynder(final double[] x, final double[] y, int n) {
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
