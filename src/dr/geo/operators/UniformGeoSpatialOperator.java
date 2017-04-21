/*
 * UniformGeoSpatialOperator.java
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

package dr.geo.operators;

import dr.geo.AbstractPolygon2D;
import dr.geo.Polygon2D;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.awt.geom.Point2D;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class UniformGeoSpatialOperator extends SimpleMCMCOperator {

    public UniformGeoSpatialOperator(Parameter parameter, double weight, List<AbstractPolygon2D> polygonList) {
        this.parameter = parameter;
        this.polygonList = polygonList;
        setWeight(weight);
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        double[] mass = null;
        int whichRegion = 0;
        if (polygonList.size() > 1) {
            mass = new double[polygonList.size()];
            for (int i = 0; i < polygonList.size(); ++i) {
                mass[i] = polygonList.get(i).calculateArea();
            }
            whichRegion = MathUtils.randomChoicePDF(mass);
        }

        totalOps++;
        if (whichRegion == 0) {
            currentSum++; // For debugging prior probability of first region
        }

        AbstractPolygon2D polygon = polygonList.get(whichRegion);

        double[][] minMax = polygon.getXYMinMax();

        int attempts = 0;
        Point2D pt;
        do {
            pt = new Point2D.Double(
                    (MathUtils.nextDouble() * (minMax[0][1] - minMax[0][0])) + minMax[0][0],
                    (MathUtils.nextDouble() * (minMax[1][1] - minMax[1][0])) + minMax[1][0]);
            attempts++;

        } while (!polygon.containsPoint2D(pt));

        if (DEBUG) {
            System.err.println("region: " + whichRegion + " attempts: " + attempts + " " + mass[0] + " " + mass[1] + "     " + (
                    (double) currentSum / (double) totalOps
            ));
        }

        parameter.setParameterValue(0, pt.getX());
        parameter.setParameterValue(1, pt.getY());

        return 0.0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "uniformGeoSpatial(" + parameter.getParameterName() + ")";
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String toString() {
        return UniformGeoSpatialOperatorParser.UNIFORM_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    //PRIVATE STUFF
    private final Parameter parameter;
    private final List<AbstractPolygon2D> polygonList;

    private long totalOps = 0;
    private long currentSum = 0;

    private static final boolean DEBUG = false;
}
