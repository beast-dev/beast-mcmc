/*
 * Polygon2DFill.java
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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.xml.XMLParseException;
import org.jdom.Element;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * 2D Polygon that uses a samplingProbability attribute as inclusion probability.
 * These polygons are typically non-overlapping and a list of such polygons currently needs its fillValues to sum to one.
 * @author Guy Baele
 */
public class Polygon2DSampling extends Polygon2D {

    public Polygon2DSampling(double[] x, double[] y, double fillValue) {
        super(x,y);
        this.fillValue = fillValue;
    }

    public Polygon2DSampling(Element e, double fillValue) {
        super(e);

        this.fillValue = fillValue;
        //no need to recompute this all the time
        this.logFillValue = Math.log(fillValue);
        //this.hasFillValue = true;
    }

    public double getProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (contains) {
            return this.fillValue;
        } else {
            return 0.0;
        }
    }

    public double getLogProbability(Point2D Point2D, boolean outside) {
        boolean contains = containsPoint2D(Point2D);
        if (contains) {
            return this.logFillValue;
        } else {
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public void setFillValue(double value) {
        fillValue = value;
        logFillValue = Math.log(value);
        //hasFillValue = true;
    }

    @Override
    public double getFillValue() {
        return fillValue;
    }

    @Override
    public double getLogFillValue() {
        return logFillValue;
    }

    @Override
    public boolean hasFillValue() {
        return true;
    }

    private double fillValue;
    private double logFillValue;
    //private boolean hasFillValue;

    public static void main(String[] args) {

        List<AbstractPolygon2D> polygons = Polygon2D.readKMLFile(args[0]);
        System.out.println("Reading from file: " + args[0] + "\n");

        List<GeoSpatialDistribution> geoSpatialDistributions = new ArrayList<GeoSpatialDistribution>();
        //check if all the polygons have either a fillValue or do not have this fillValue
        boolean overAllFill = false;
        AbstractPolygon2D first = polygons.get(0);
        System.out.println("Polygon 0: " + first.getFillValue());
        boolean equalFills = first.hasFillValue();
        if (polygons.size() > 1) {
            for (int i = 1; i < polygons.size(); i++) {
                if (!(equalFills == polygons.get(i).hasFillValue())) {
                    throw new RuntimeException("Inconsistent fillValue attributes provided.");
                } else {
                    System.out.println("Polygon " + i + ": " + polygons.get(i).getFillValue());
                }
            }
        }
        if (equalFills) {
            double sum = 0.0;
            for (AbstractPolygon2D region : polygons) {
                sum += region.getFillValue();
            }
            if (Math.abs(sum - 1.0) > 1E-12) {
                throw new RuntimeException("Fill values in " + args[0] + " do not sum to 1 : " + sum);
            } else {
                System.out.println("Provided sampling probabilities sum to 1.0");
            }
        }
        overAllFill = equalFills;
        for (AbstractPolygon2D region : polygons) {
            geoSpatialDistributions.add(new GeoSpatialDistribution("", region, true));
        }

        System.out.println("\nTesting sampling probabilities:");

        MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
                new MultiRegionGeoSpatialDistribution("", geoSpatialDistributions, false, overAllFill));
        double[] coordinates = new double[2];
        coordinates[0] = 12.68379879;
        coordinates[1] = 108.12982178;
        Parameter parameter = new Parameter.Default(coordinates);
        likelihood.addData(parameter);

        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        parameter.setParameterValue(0,13.68379879);
        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        parameter.setParameterValue(0,12.68);
        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        parameter.setParameterValue(0, 13.03387928);
        parameter.setParameterValue(1, 108.51604462);
        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        parameter.setParameterValue(0, 13.01971245);
        parameter.setParameterValue(1, 108.56764221);
        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        parameter.setParameterValue(0, 13.018);
        System.out.println("logL = " + likelihood.getLogLikelihood());
        System.out.println("L = " + Math.exp(likelihood.getLogLikelihood()));

        //from 8 to 14 and from 104 to 110
        /*System.out.println("\n\n\nplot(0,0,xlim=c(9,14),ylim=c(105,110))");
        Random rand = new Random();
        for (int i = 0; i < 100000; i++) {
            parameter.setParameterValue(0, rand.nextDouble()*6+9);
            parameter.setParameterValue(1, rand.nextDouble()*6+105);
            if (!Double.isInfinite(likelihood.getLogLikelihood())) {
                double color = Math.exp(likelihood.getLogLikelihood());
                //System.out.println(color);
                String col = "black";
                if (Math.abs(color - 0.106508875739645) < 1E-5) {
                    col = "red";
                } else if (Math.abs(color - 0.029585798816568) < 1E-5) {
                    col = "blue";
                } else if (Math.abs(color - 0.0769230769230769) < 1E-5) {
                    col = "green";
                } else if (Math.abs(color - 0.0355029585798817) < 1E-5) {
                    col = "grey";
                } else if (Math.abs(color - 0.124260355029586) < 1E-5) {
                    col = "purple";
                } else if (Math.abs(color - 0.112426035502959) < 1E-5) {
                    col = "orange";
                } else if (Math.abs(color - 0.136094674556213) < 1E-5) {
                    col = "yellow";
                } else if (Math.abs(color - 0.0946745562130177) < 1E-5) {
                    col = "burlywood";
                } else if (Math.abs(color - 0.284023668639053) < 1E-5) {
                    col = "deepskyblue";
                }
                System.out.println("points(" + parameter.getParameterValue(0) + "," + parameter.getParameterValue(1) + ",col=\"" + col + "\",pch=20,bg=\""+ col +"\")");
            }
        }*/

    }

}
