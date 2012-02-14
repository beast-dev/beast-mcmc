package dr.inference.distribution;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.InverseGaussianDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.InverseGaussianDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/*
 * InverseGaussianDistributionModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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


/**
 * @author Wai Lok Sibon Li
 * @version $Id: InverseGaussianDistributionModel.java,v 1.8 2009/03/30 20:25:59 rambaut Exp $
 */

public class InverseGaussianDistributionModel extends AbstractModel implements ParametricDistributionModel {

    /**
     * @param meanParameter  the mean, mu
     * @param igParameter   either the standard deviation parameter, sigma or the shape parameter, lamba
     * @param offset         offset of the distribution
     * @param useShape         whether shape or stdev is used
     */
    public InverseGaussianDistributionModel(Parameter meanParameter, Parameter igParameter, double offset, boolean useShape) {

        super(InverseGaussianDistributionModelParser.INVERSEGAUSSIAN_DISTRIBUTION_MODEL);

        if(useShape) {
            this.shapeParameter = igParameter;
            this.stdevParameter = null;
            addVariable(shapeParameter);
            this.shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
        else {
            this.stdevParameter = igParameter;
            this.shapeParameter = null;
            addVariable(stdevParameter);
            this.stdevParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        this.meanParameter = meanParameter;
        addVariable(meanParameter);
        this.offset = offset;
        this.meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

    }

    public final double getS() {
        if(stdevParameter==null) {
            return Math.sqrt(InverseGaussianDistribution.variance(getM(), getShape()));
        }
        return stdevParameter.getParameterValue(0);
    }

    public final void setS(double S) {
        if(stdevParameter==null) {
            throw new RuntimeException("Standard deviation parameter is not being used");
        }
        else {
            stdevParameter.setParameterValue(0, S);
        }
    }

    public final Parameter getSParameter() {
        if(stdevParameter==null) {
            throw new RuntimeException("Standard deviation parameter is not being used");
        }
        return stdevParameter;
    }

    public final double getShape() {
        if(shapeParameter == null) {
            double shape = (getM() * getM() * getM()) / (getS() * getS());
            return shape;
        }
        return shapeParameter.getParameterValue(0);
    }

    public final void setShape(double shape) {
        if(shapeParameter==null) {
            throw new RuntimeException("Shape parameter is not being used");

        }
        else {
            shapeParameter.setParameterValue(0, shape);
        }
    }

    public final Parameter getShapeParameter() {
        if(shapeParameter==null) {
            throw new RuntimeException("Shape parameter is not being used");
        }
        return shapeParameter;
    }

    /* Unused method */
    //private double getStDev() {
    //return Math.sqrt(InverseGaussianDistribution.variance(getM(), getShape()));//Math.sqrt((getM()*getM()*getM())/getShape());
    //}

    /**
     * @return the mean
     */
    public final double getM() {
        return meanParameter.getParameterValue(0);
    }

    public final void setM(double M) {
        meanParameter.setParameterValue(0, M);
        //double shape = (getM() * getM() * getM()) / (getS() * getS());
        //setShape(shape);
    }

    public final Parameter getMParameter() {
        return meanParameter;
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return InverseGaussianDistribution.pdf(x - offset, getM(), getShape());
    }

    public double logPdf(double x) {
        if (x - offset <= 0.0) return Double.NEGATIVE_INFINITY;
        return InverseGaussianDistribution.logPdf(x - offset, getM(), getShape());
    }

    public double cdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return InverseGaussianDistribution.cdf(x - offset, getM(), getShape());
    }

    public double quantile(double y) {
        return InverseGaussianDistribution.quantile(y, getM(), getShape()) + offset;
    }

    /**
     * @return the mean of the distribution
     */
    public double mean() {
        //return InverseGaussianDistribution.mean(getM(), getShape()) + offset;
        return getM() + offset;
    }

    /**
     * @return the variance of the distribution.
     */
    public double variance() {
        //return InverseGaussianDistribution.variance(getM(), getShape());
        return getS() * getS();
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            System.out.println("just checking if this ever gets used anyways... probably have to change the getLowerBound in LogNormalDistributionModel if it does");
            return pdf(x);
        }

        public final double getLowerBound() {
            return 0.0;
            //return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented!");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Parameter meanParameter;
    private final Parameter stdevParameter;
    private final Parameter shapeParameter;
    private final double offset;

}