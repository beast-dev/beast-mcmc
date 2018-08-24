/*
 * DistributionLikelihood.java
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

import dr.math.distributions.Distribution;
import dr.math.matrixAlgebra.Vector;
import dr.util.Attribute;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to the given parametric distribution.
 *
 * @author Alexei Drummond
 * @version $Id: DistributionLikelihood.java,v 1.11 2005/05/25 09:35:28 rambaut Exp $
 */

public class DistributionLikelihood extends AbstractDistributionLikelihood {

    public final static boolean DEBUG = false;

    public static final String DISTRIBUTION_LIKELIHOOD = "distributionLikelihood";

    private int from = -1;
    private int to = Integer.MAX_VALUE;
    private final boolean evaluateEarly;

    public DistributionLikelihood(Distribution distribution) {
        this(distribution, 0.0, false, 1.0);
    }

    public DistributionLikelihood(Distribution distribution, double offset) {
        this(distribution, offset, offset > 0.0, 1.0);
    }

    public DistributionLikelihood(Distribution distribution, double offset, double scale){
        this(distribution, offset, offset>0.0, scale);
    }

    public DistributionLikelihood(Distribution distribution, boolean evaluateEarly) {
        this(distribution, 0.0, evaluateEarly, 1.0);
    }

    public DistributionLikelihood(Distribution distribution, double offset, boolean evaluateEarly, double scale) {
        super(null);
        this.distribution = distribution;
        this.offset = offset;
        this.evaluateEarly = evaluateEarly;
        this.scale=scale;
    }

    public DistributionLikelihood(ParametricDistributionModel distributionModel) {
        super(distributionModel);
        this.distribution = distributionModel;
        this.offset = 0.0;
        this.evaluateEarly = false;
        this.scale=1.0;
    }

    public Distribution getDistribution() {
        return distribution;
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

        if (DEBUG) {
            System.err.println("Calling DistributionLikelihood.calculateLogLikelihood()");
            System.err.println(distribution.toString());
            System.err.println(dataList.toString() + "\n");
        }

        double logL = 0.0;
        int count = 0;

        for( Attribute<double[]> data : dataList ) {

            // Using this in the loop is incredibly wasteful, especially in the loop condition to get the length
            final double[] attributeValue = data.getAttributeValue();

            if (DEBUG) {
                System.err.println("\t" + new Vector(attributeValue));
            }

            for (int j = Math.max(0, from); j < Math.min(attributeValue.length, to); j++) {

                final double value = attributeValue[j] - offset;

                if (offset > 0.0 && value < 0.0) {
                    // fixes a problem with the offset on exponential distributions not
                    // actually bounding the distribution. This only performs this check
                    // if a non-zero offset is actually given otherwise it assumes the
                    // parameter is either legitimately allowed to go negative or is bounded
                    // at zero anyway.
                    return Double.NEGATIVE_INFINITY;
                }
                logL += getLogPDF(value, count);
                count += 1;
            }

        }
        return logL;
    }

    protected double getLogPDF(double value, int i){
        return distribution.logPdf(value/scale) - Math.log(scale);
    }

    @Override
    public boolean evaluateEarly() {
        return evaluateEarly;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    @Override
    public String prettyName() {
        String s = distribution.getClass().getName();
        String[] parts = s.split("\\.");
        s = parts[parts.length - 1];
        if( s.endsWith("Distribution") ) {
            s = s.substring(0, s.length() - "Distribution".length());
        }
        s = s + '(';
        for( Attribute<double[]> data : dataList ) {
            String name = data.getAttributeName();
            if( name == null ) {
                name = "?";
            }
                s = s + name + ',';
        }
        s = s.substring(0,s.length()-1) + ')';

        return s;
    }

    protected Distribution distribution;
    private final double offset;
    private final double scale;
}

