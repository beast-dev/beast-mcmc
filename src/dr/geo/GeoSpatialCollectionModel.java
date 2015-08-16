/*
 * GeoSpatialCollectionModel.java
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

package dr.geo;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

/**
 * @author Marc A. Suchard
 *         <p/>
 *         Provides a GeoSpatialDistribution over multiple points in multiple polygon.
 *         Uses AbstractModelLikelihood to cache 'contains' to reduce recalculations
 *         when only a single point is updated
 */

public class GeoSpatialCollectionModel extends AbstractModelLikelihood {

    public GeoSpatialCollectionModel(String name, Parameter points,
                                     List<GeoSpatialDistribution> geoSpatialDistributions,
                                     boolean isIntersection) {

        super(name);
        this.points = points;
        this.geoSpatialDistributions = geoSpatialDistributions;

        dim = points.getDimension() / GeoSpatialDistribution.dimPoint;
        cachedPointLogLikelihood = new double[dim];
        storedCachedPointLogLikelihood = new double[dim];
        validPointLogLikelihood = new boolean[dim];
        storedValidPointLogLikelihood = new boolean[dim];
        likelihoodKnown = false;

        addVariable(points);

        this.isIntersection = isIntersection;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No submodels; do nothing
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Mark appropriate dim as invalid
        validPointLogLikelihood[index / GeoSpatialDistribution.dimPoint] = false;
        likelihoodKnown = false;
    }

    protected void storeState() {

        System.arraycopy(cachedPointLogLikelihood, 0, storedCachedPointLogLikelihood, 0, dim);
        System.arraycopy(validPointLogLikelihood, 0, storedValidPointLogLikelihood, 0, dim);

        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {

        double[] tmp1 = cachedPointLogLikelihood;
        cachedPointLogLikelihood = storedCachedPointLogLikelihood;
        storedCachedPointLogLikelihood = tmp1;

        boolean[] tmp2 = validPointLogLikelihood;
        validPointLogLikelihood = storedValidPointLogLikelihood;
        storedValidPointLogLikelihood = tmp2;

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected void acceptState() {

    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {

        if (likelihoodKnown)
            return logLikelihood;

        logLikelihood = 0.0;
        final double[] point = new double[GeoSpatialDistribution.dimPoint];

        for (int i = 0; i < dim; i++) {
            if (!validPointLogLikelihood[i]) {
                final int offset = i * GeoSpatialDistribution.dimPoint;
                for (int j = 0; j < GeoSpatialDistribution.dimPoint; j++)
                    point[j] = points.getParameterValue(offset + j);

                double pointLogLikelihood = 0;
                for (GeoSpatialDistribution distribution : geoSpatialDistributions) {
                    //if we consider the union of polygons and the point must be inside, than it is good enough that the point is in one polygon
                    //so we look for a polygon that does not yield -inf
                    if (!isIntersection && !distribution.getOutside()) {
                        final double logPdf = distribution.logPdf(point);
                        if (logPdf != Double.NEGATIVE_INFINITY) {
                            pointLogLikelihood = logPdf;
                            break;
                        } else {
                            pointLogLikelihood = logPdf;
                        }
                    } else {
                        // Below is for intersections (or unions of complements)
                        pointLogLikelihood += distribution.logPdf(point);
                        if (pointLogLikelihood == Double.NEGATIVE_INFINITY)
                            break; // No need to finish
                    }
                }
                cachedPointLogLikelihood[i] = pointLogLikelihood;
                validPointLogLikelihood[i] = true;
            }
            logLikelihood += cachedPointLogLikelihood[i];
            if (logLikelihood == Double.NEGATIVE_INFINITY)
                break; // No need to finish
        }
        likelihoodKnown = true;
        return logLikelihood;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        for (int i = 0; i < dim; i++)
            validPointLogLikelihood[i] = false;
    }

    public Parameter getParameter() {
        return points;
    }

    private Parameter points;
    private List<GeoSpatialDistribution> geoSpatialDistributions;
    private int dim;

    private double[] cachedPointLogLikelihood;
    private double[] storedCachedPointLogLikelihood;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean[] validPointLogLikelihood;
    private boolean[] storedValidPointLogLikelihood;

    private final boolean isIntersection;
}
