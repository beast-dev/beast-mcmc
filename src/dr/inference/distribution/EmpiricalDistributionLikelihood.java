/*
 * EmpiricalDistributionLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.GradientProvider;
import dr.util.Attribute;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to a distribution generated empirically from some data.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id:$
 */

public abstract class EmpiricalDistributionLikelihood extends AbstractDistributionLikelihood
        implements GradientProvider {

    public static final String EMPIRICAL_DISTRIBUTION_LIKELIHOOD = "empiricalDistributionLikelihood";

    private static final double MIN_DENSITY_PROPORTION = 1E-2;
    private static final boolean DEBUG = false;

    private int from = -1;
    private int to = Integer.MAX_VALUE;
    private double offset = 0;
    private double lower = Double.NEGATIVE_INFINITY;
    private double upper = Double.POSITIVE_INFINITY;

    EmpiricalDistributionLikelihood(String fileName, boolean inverse, boolean byColumn) {
        super(null);
        this.fileName = fileName;

        boolean densityInLogSpace = false; // TODO Read from file

        if (byColumn)
            readFileByColumn(fileName, densityInLogSpace);
        else
            readFileByRow(fileName, densityInLogSpace);

        this.inverse = inverse;    
    }

    EmpiricalDistributionLikelihood(List<EmpiricalDistributionData> densityList, boolean inverse) {
        super(null);
        this.densityList = densityList;
        this.inverse = inverse;
    }

    public void setBounds(double lower, double upper) {
        this.lower = lower;
        this.upper = upper;
    }

    class ComparablePoint2D extends Point2D.Double implements Comparable<ComparablePoint2D> {

        ComparablePoint2D(double x, double y) {
            super(x,y);
        }

        public int compareTo(ComparablePoint2D pt0) {
            return java.lang.Double.compare(getX(), pt0.getX());
        }
    }

    private void readFileByRow(String fileName, boolean densityInLogSpace) {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));

            List<ComparablePoint2D> ptList = new ArrayList<>();

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.charAt(0) != '#') {
                    StringTokenizer st = new StringTokenizer(line);
                    try {
                        double value = Double.valueOf(st.nextToken());
                        double density = Double.valueOf(st.nextToken());
                        ptList.add(new ComparablePoint2D(value,density));
                    } catch (Exception e) {
                        System.err.println("Error parsing line: '"+line+"' in "+fileName);
                        System.exit(-1);
                    }
                }
            }
            Collections.sort(ptList);
            // Prune off begining and ending zeros

            while( (ptList.get(0).getY() == 0) &&
                   (ptList.get(1).getY() == 0) )
                ptList.remove(0);

            while( (ptList.get(ptList.size()-1).getY() == 0) &&
                   (ptList.get(ptList.size()-2).getY() == 0) )
                ptList.remove(ptList.size()-1);
            
            // Find min density
            double minDensity = Double.POSITIVE_INFINITY;
            for(ComparablePoint2D pt : ptList) {
                if (pt.getY() > 0.0 && pt.getY() < minDensity)
                    minDensity = pt.getY();
            }
            // Set zeros in the middle to 1/100th of minDensity
            for(ComparablePoint2D pt : ptList) {
                if (pt.getY() == 0)
                    pt.y = minDensity * MIN_DENSITY_PROPORTION;
            }
            EmpiricalDistributionData data = new EmpiricalDistributionData(
                    new double[ptList.size()], new double[ptList.size()], densityInLogSpace);
            double total = 0.0;
            for(int i=0; i<ptList.size(); i++) {
                ComparablePoint2D pt = ptList.get(i);
                data.values[i] = pt.getX();
                data.density[i] = pt.getY();
                total += pt.getY();
            }
            for (int i = 0; i < data.density.length; ++i) {
                data.density[i] /= total;
            }
            reader.close();

            if (DEBUG) {
                System.err.println("EDL File : " + fileName);
                System.err.println("Min value: " + data.values[0]);
                System.err.println("Max value: " + data.values[data.values.length-1]);
            }

            this.densityList = new ArrayList<>();
            this.densityList.add(data);

        } catch (FileNotFoundException e) {
            System.err.println("File not found: "+fileName);
            System.exit(-1);
        } catch (IOException e) {
            System.err.println("IO exception reading: "+fileName);
            System.exit(-1);
        }
    }

    private void readFileByColumn(String fileName, boolean densityInLogSpace) {

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));

            String line1 = reader.readLine();
            StringTokenizer st = new StringTokenizer(line1," ");
            double[] values = new double[st.countTokens()];
            for(int i=0; i<values.length; i++)
                values[i] = Double.valueOf(st.nextToken());
            String line2 = reader.readLine();
            st = new StringTokenizer(line2," ");
            double[] density = new double[st.countTokens()];
            for(int i=0; i<density.length; i++)
                density[i] = Double.valueOf(st.nextToken());

            densityList = new ArrayList<>();
            densityList.add(new EmpiricalDistributionData(values, density, densityInLogSpace));

            reader.close();

        } catch (FileNotFoundException e) {
            System.err.println("File not found: "+fileName);
            System.exit(-1);
        } catch (IOException e) {
            System.err.println("IO exception reading: "+fileName);
            System.exit(-1);
        }

    }

    public void setOffset(double offset) {
        this.offset = offset;
    }


    public void setRange(int from, int to) {
        this.from = from;
        this.to = to;
    }
  
    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;

        for (Attribute<double[]> data : dataList) {  
            // see comment in DistributionLikelihood
            final double[] attributeValue = data.getAttributeValue();

            int index = 0;
            for (int j = Math.max(0, from); j < Math.min(attributeValue.length, to); j++) {

                double value = attributeValue[j] + offset;
                if (value > lower && value < upper) {
                    logL += logPDF(value, densityList.get(index));
                } else {
                    return Double.NEGATIVE_INFINITY;
                }

                if (DEBUG) {
                    if (Double.isInfinite(logL)) {
                        System.err.println("Infinite log density in "+
                                (getId() != null ? getId() : fileName)                        
                        );
                        System.err.println("Evaluated at "+value);
                        EmpiricalDistributionData d = densityList.get(index);
                        System.err.println("Min: " + d.values[0] + " Max: " + d.values[d.values.length - 1]);
                    }
                }

                if (densityList.size() > 1) {
                    ++index;
                }
            }
        }
        return logL;
    }

    public int getDistributionDimension() {
        return densityList.size();
    }

    public int getDimension() {

        int size = 0;
        for (Attribute<double[]> data : dataList) {
            size += data.getAttributeValue().length;
        }

        return size;
    }

    public double[] getGradientLogDensity(double[] point) {

        assert densityList.size() == 1 || densityList.size() == point.length;

        double[] gradient = new double[point.length];

        int index = 0;
        for (int i = 0; i < point.length; ++i) {
            double x = point[i] + offset;

            gradient[i] = (x > lower && x < upper) ? gradientLogPdf(x, densityList.get(index)) : 0.0;

            if (densityList.size() > 1) {
                ++index;
            }
        }

        return gradient;
    }

    public double[] getGradientLogDensity(Object x) {
        return getGradientLogDensity(GradientProvider.toDoubleArray(x));
    }

    abstract protected double logPDF(double value, EmpiricalDistributionData data);

    abstract protected double gradientLogPdf(double value, EmpiricalDistributionData data);

    private List<EmpiricalDistributionData> densityList;
//    protected Data data;

//    protected  double[] values;
//    protected  double[] density;

    protected boolean inverse;

    protected String fileName;
}

