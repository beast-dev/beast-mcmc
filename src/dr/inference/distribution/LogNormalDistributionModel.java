/*
 * LogNormalDistributionModel.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.NormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for log-normally distributed data.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: LogNormalDistributionModel.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */

public class LogNormalDistributionModel extends AbstractModel implements ParametricDistributionModel {

    public enum Parameterization {
        MU_SIGMA,
        MU_PRECISION,
        MEAN_STDEV
    }

    /**
     * Constructor.
     * This is the old constructor left for backwards compatibility
     */
    public LogNormalDistributionModel(Parameter meanParameter, Parameter stdevParameter, double offset, boolean parametersInRealSpace) {

        super(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);

        this.offset = offset;

        if (parametersInRealSpace) {
            this.meanParameter = meanParameter;
            this.muParameter = null;
            addVariable(this.meanParameter);
            this.meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

            this.stdevParameter = stdevParameter;
            this.sigmaParameter = null;
            addVariable(this.stdevParameter);
            this.stdevParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

            this.parameterization = Parameterization.MEAN_STDEV;
        } else {
            this.muParameter = meanParameter;
            this.meanParameter = null;
            addVariable(this.muParameter);
            this.muParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

            this.sigmaParameter = stdevParameter;
            this.stdevParameter = null;
            addVariable(this.sigmaParameter);
            this.sigmaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

            this.parameterization = Parameterization.MU_SIGMA;
        }

        precisionParameter = null;
    }

    /**
     * Constructor that allows different parameterizations. Generally we would expect either mu and sigma or mean and
     * stdev to be provided (the others should be null). But the only restriction is one of mean or mu and one of stdev
     * and sigma.
     * @param parameterization
     * @param parameter1
     * @param parameter2
     * @param offset
     */
    public LogNormalDistributionModel(Parameterization parameterization, Parameter parameter1, Parameter parameter2, double offset) {

        super(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL);

        switch (parameterization) {
            case MU_SIGMA:
                this.muParameter = parameter1;
                this.sigmaParameter = parameter2;
                this.meanParameter = null;
                this.stdevParameter = null;
                this.precisionParameter = null;
                this.parameterization = Parameterization.MU_SIGMA;
                break;
            case MU_PRECISION:
                this.muParameter = parameter1;
                this.precisionParameter = parameter2;
                this.meanParameter = null;
                this.stdevParameter = null;
                this.sigmaParameter = null;
                this.parameterization = Parameterization.MU_PRECISION;
                break;
            case MEAN_STDEV:
                this.meanParameter = parameter1;
                this.stdevParameter = parameter2;
                this.muParameter = null;
                this.sigmaParameter = null;
                this.precisionParameter = null;
                this.parameterization = Parameterization.MEAN_STDEV;
                break;
            default:
                throw new IllegalArgumentException("Unknow parameterization type");
        }


        this.offset = offset;

        if (this.muParameter != null) {
            addVariable(this.muParameter);
            this.muParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        }

        if (meanParameter != null) {
            addVariable(this.meanParameter);
            this.meanParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        if (sigmaParameter != null) {
            addVariable(this.sigmaParameter);
            this.sigmaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        if (stdevParameter != null) {
            addVariable(this.stdevParameter);
            this.stdevParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        if (precisionParameter != null) {
            addVariable(this.precisionParameter);
            this.precisionParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
    }

    public Parameter getMeanParameter() {
        return meanParameter;
    }

    public Parameter getStdevParameter() {
        return stdevParameter;
    }

    public Parameter getMuParameter() {
        return muParameter;
    }

    public Parameter getSigmaParameter() {
        return sigmaParameter;
    }

    public Parameter getPrecisionParameter() {
        return precisionParameter;
    }

    public Parameterization getParameterization() {
        return parameterization;
    }

    public final double getMu() {
        if (muParameter != null) {
            return muParameter.getValue(0);
        } else {
            return calculateMu(meanParameter.getValue(0), stdevParameter.getValue(0));
        }
    }

    public final double getMean() {
        if (meanParameter != null) {
            return meanParameter.getValue(0);
        } else {
            return calculateMean(muParameter.getValue(0), getSigma());
        }
    }

    public final double getSigma() {
        if (sigmaParameter != null) {
            return sigmaParameter.getValue(0);
        } else if (precisionParameter != null) {
            return Math.sqrt(1.0 / precisionParameter.getValue(0));
        } else{
            return calculateSigma(meanParameter.getValue(0), stdevParameter.getValue(0));
        }

    }

    public final double getStdev() {
        if (stdevParameter != null) {
            return stdevParameter.getValue(0);
        } else if (precisionParameter != null) {
            double sigma = Math.sqrt(1.0 / precisionParameter.getValue(0));
            return calculateStdev(muParameter.getValue(0), sigma);
        } else {
            return calculateStdev(muParameter.getValue(0), sigmaParameter.getValue(0));
        }
    }

    public final double getPrecision() {
        if (precisionParameter != null) {
            return precisionParameter.getValue(0);
        } else {
            double sigma = getSigma();
            return 1.0 / (sigma * sigma);
        }
    }

    private double calculateMu(double mean, double stdev) {
        return Math.log(mean/Math.sqrt(1 + (stdev * stdev) / (mean * mean)));
    }

    private double calculateSigma(double mean, double stdev) {
        return Math.sqrt(Math.log(1 + (stdev * stdev) / (mean * mean)));
    }

    private double calculateMean(double mu, double sigma) {
        return Math.exp(mu + 0.5 * sigma * sigma);
    }

    private double calculateStdev(double mu, double sigma) {
        return Math.sqrt( (Math.exp(sigma * sigma)-1) * Math.exp(2.0 * mu + (sigma * sigma)) );
    }

//    public final double getS() {
//        //System.out.println(isStdevInRealSpace+"\t" + isMeanInRealSpace + "\t" + Math.sqrt(Math.log(1 + scaleParameter.getParameterValue(0)/Math.pow(meanParameter.getParameterValue(0), 2))) + "\t" + scaleParameter.getParameterValue(0));
//        if(isStdevInRealSpace) {
//
//            if(isMeanInRealSpace) {
//                return Math.sqrt(Math.log(1 + scaleParameter.getParameterValue(0)/Math.pow(meanParameter.getParameterValue(0), 2)));
//            }
//            else {
//                throw new RuntimeException("S can not be computed with M and stdev");
//            }
//        }
//        return scaleParameter.getParameterValue(0);
//    }
//
//    public final void setS(double S) {
//        scaleParameter.setParameterValue(0, S);
//    }
//
//    public final Parameter getSParameter() {
//        return scaleParameter;
//    }
//
//    /* StDev in this class is actually incorrectly named the S parameter */
//    private double getStDev() {
//        return usesStDev ? getS() : Math.sqrt(1.0 / getS());
//    }
//
//    /**
//     * @return the mean (always in log space)
//     */
//    public final double getM() {
//        if (isMeanInRealSpace) {
//            double stDev = getStDev();
//            return Math.log(meanParameter.getParameterValue(0)) - (0.5 * stDev * stDev);
//        } else {
//            return meanParameter.getParameterValue(0);
//
//        }
//    }
//
//    public final void setM(double M) {
//        if (isMeanInRealSpace) {
//            double stDev = getStDev();
//            meanParameter.setParameterValue(0, Math.exp(M + (0.5 * stDev * stDev)));
//        } else {
//            meanParameter.setParameterValue(0, M);
//        }
//    }
//
//    public final Parameter getMeanParameter() {
//        return meanParameter;
//    }
//
//     public Parameter getPrecisionParameter() {
//        if (!usesStDev)
//            return scaleParameter;
//        return null;
//    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return NormalDistribution.pdf(Math.log(x - offset), getMu(), getSigma()) / (x - offset);
    }


    public double logPdf(double x) {
        if (x - offset <= 0.0) return Double.NEGATIVE_INFINITY;
        return NormalDistribution.logPdf(Math.log(x - offset), getMu(), getSigma()) - Math.log(x - offset);
    }

    public double cdf(double x) {
        if (x - offset <= 0.0) return 0.0;
        return NormalDistribution.cdf(Math.log(x - offset), getMu(), getSigma());
    }

    public double quantile(double y) {
        return Math.exp(NormalDistribution.quantile(y, getMu(), getSigma())) + offset;
    }

    /**
     * @return the mean of the distribution
     */
    public double mean() {
        return getMean();
    }

    /**
     * @return the variance of the log normal distribution.  Not really the variance of the lognormal but the S^2
     * parameter
     */
    public double variance() {
        return getStdev() * getStdev();
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        return pdfFunction;
    }

    private final UnivariateFunction pdfFunction = new UnivariateFunction() {
        public final double evaluate(double x) {
            return pdf(Math.log(x));
        }

        public final double getLowerBound() {
            return Double.NEGATIVE_INFINITY;
        }

        public final double getUpperBound() {
            return Double.POSITIVE_INFINITY;
        }
    };

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return meanParameter;
    }

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

    // Can be parameterized as the mean and stdev of the lognormal or the mean and stdev of the underlying normal (mu & sigma)
    private final Parameter meanParameter;
    private final Parameter stdevParameter;
    private final Parameter muParameter;
    private final Parameter sigmaParameter;
    private final Parameter precisionParameter;
    private final double offset;

    private Parameterization parameterization;

    public static void main(String[] argv) {
        Parameter meanParameter = new Parameter.Default(1.0);
        Parameter stdevParameter = new Parameter.Default(5.0);
        Parameter muParameter = new Parameter.Default(-1.629048);
        Parameter sigmaParameter = new Parameter.Default(1.80502);

        LogNormalDistributionModel ln1 = new LogNormalDistributionModel(Parameterization.MEAN_STDEV, meanParameter, stdevParameter, 0);
        System.out.println("Lognormal mean = 1.0, stdev = 5.0");
        System.out.println("  mu = " + ln1.getMu() + " (correct = -1.629048)");
        System.out.println("  sigma = " + ln1.getSigma() + " (correct = 1.80502)");
        System.out.println("  quantile(2.5) = " + ln1.quantile(0.025) + " (correct = 0.005702663)");
        System.out.println("  quantile(97.5) = " + ln1.quantile(0.975) + " (correct = 6.744487892)");

        LogNormalDistributionModel ln2 = new LogNormalDistributionModel(Parameterization.MU_SIGMA, muParameter, sigmaParameter, 0);
        System.out.println("Lognormal mu = -1.629048, sigma = 1.80502");
        System.out.println("  mean = " + ln2.getMean() + " (correct = 1.0)");
        System.out.println("  sigma = " + ln2.getStdev() + " (correct = 5.0)");
        System.out.println("  quantile(2.5) = " + ln2.quantile(0.025) + " (correct = 0.005702663)");
        System.out.println("  quantile(97.5) = " + ln2.quantile(0.975) + " (correct = 6.744487892)");

        meanParameter = new Parameter.Default(0.001);
        stdevParameter = new Parameter.Default(0.0005);
        LogNormalDistributionModel ln3 = new LogNormalDistributionModel(Parameterization.MEAN_STDEV, meanParameter, stdevParameter, 0);
        System.out.println("Lognormal mean = 0.001, stdev = 0.0005");
        System.out.println("  mu = " + ln3.getMu());
        System.out.println("  sigma = " + ln3.getSigma());
        for (int i = 1; i <= 12; i++) {
            double y = ((double)i) / 13.0;
            System.out.println(i + "\t" + y + "\t" + ln3.quantile(y));
        }


    }
}
