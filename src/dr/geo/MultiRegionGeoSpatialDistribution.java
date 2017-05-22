/*
 * MultiRegionGeoSpatialDistribution.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.geo;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class MultiRegionGeoSpatialDistribution extends GeoSpatialDistribution {

    public MultiRegionGeoSpatialDistribution(String label) {
        super(label);
        regions = new ArrayList<GeoSpatialDistribution>();
        union = false;
        fillValue = false;
    }

    public MultiRegionGeoSpatialDistribution(String label, List<GeoSpatialDistribution> regions, boolean union, boolean fillValue) {
        super(label);
        this.regions = regions;
        this.union = union;
        this.fillValue = fillValue;
    }

    public double logPdf(double[] x) {

        if (fillValue) {
            //System.err.println("fillValue set to true");
            //int test = 0;
            for (GeoSpatialDistribution region : regions) {
                //System.err.println(test);
                //test++;
                if (!Double.isInfinite(region.logPdf(x))) {
                    return region.logPdf(x);
                }
            }
            return Double.NEGATIVE_INFINITY;
        }

        if (union) {
            for (GeoSpatialDistribution region : regions) {
                if (region.logPdf(x) == 0.0) { // matches
                    return 0.0;
                }
            }
            return Double.NEGATIVE_INFINITY;
        } // else is intersection

        for (GeoSpatialDistribution region : regions) {
            if (region.logPdf(x) == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        return 0.0;
    }

    public double[][] getScaleMatrix() {
        return null;
    }

    public double[] getMean() {
        return null;
    }

    public String getType() {
        return TYPE;
    }

    public String getLabel() {
        return label;
    }

    public boolean getUnion() {
        return union;
    }

    public List<GeoSpatialDistribution> getRegions() {
        return regions;
    }

    private final List<GeoSpatialDistribution> regions;
    private final boolean union;
    private final boolean fillValue;

}
