/*
 * NormalPeriodPriorDistribution.java
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

package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.NormalGammaDistribution;
import dr.math.functionEval.GammaFunction;
import dr.xml.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.TDistributionImpl;

import java.util.ArrayList;
import java.util.Arrays;

/**
 The assumption here is that the periods are drawn from a normal distribution with unknown mean and variance.
 The hyperprior is the conjugate, normal-gamma distribution.

 @author Matthew Hall
 */
public class NormalPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String NORMAL = "normalPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";
    public static final String MU = "mu";
    public static final String LAMBDA = "lambda";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    private NormalGammaDistribution hyperprior;

    private Parameter posteriorMean;
    private Parameter posteriorBeta;
    private Parameter posteriorExpectedPrecision;

    double normalApproximationThreshold = 30;

    private ArrayList<Double> dataValues;
    private double[] currentParameters;

    public NormalPeriodPriorDistribution(String name, boolean log, NormalGammaDistribution hyperprior){
        super(name, log);
        this.hyperprior = hyperprior;
        posteriorBeta = new Parameter.Default(1);
        posteriorMean = new Parameter.Default(1);
        posteriorExpectedPrecision = new Parameter.Default(1);
        addVariable(posteriorBeta);
        addVariable(posteriorMean);
        addVariable(posteriorExpectedPrecision);
    }

    public NormalPeriodPriorDistribution(String name, boolean log, double mu_0, double lambda_0,
                                         double alpha_0, double beta_0){
        this(name, log, new NormalGammaDistribution(mu_0, lambda_0, alpha_0, beta_0));
        reset();
    }

    public void reset(){
        dataValues = new ArrayList<Double>();
        currentParameters = hyperprior.getParameters();
        logL = 0;
    }

    // this returns the posterior predictive probability of the new value, and updates the total

    public double calculateLogPosteriorProbability(double newValue, double minValue){
        double out = calculateLogPosteriorPredictiveProbability(newValue);
        if(minValue != Double.NEGATIVE_INFINITY){
            out -= calculateLogPosteriorPredictiveCDF(minValue, true);
        }
        logL += out;
        update(newValue);
        return out;
    }

    public double calculateLogPosteriorCDF(double limit, boolean upper) {
        return calculateLogPosteriorPredictiveCDF(limit, upper);
    }


    public double calculateLogPosteriorPredictiveProbability(double value){
        double mean = currentParameters[0];
        double sd = Math.sqrt(currentParameters[3]*(currentParameters[1]+1)
                /(currentParameters[2]*currentParameters[1]));
        double scaledValue = (value - mean)/sd;
        double out;

        if(2*currentParameters[2]<=normalApproximationThreshold) {
            TDistributionImpl tDist = new TDistributionImpl(2 * currentParameters[2]);

            out = Math.log(tDist.density(scaledValue));


        } else {

            out = NormalDistribution.logPdf(scaledValue, 0, 1);

        }

        return out;
    }

    public double calculateLogPosteriorPredictiveCDF(double value, boolean upperTail){
        double mean = currentParameters[0];
        double sd = Math.sqrt(currentParameters[3]*(currentParameters[1]+1)
                /(currentParameters[2]*currentParameters[1]));
        double scaledValue = (value - mean)/sd;
        double out;

        if(2*currentParameters[2]<=normalApproximationThreshold) {
            TDistributionImpl tDist = new TDistributionImpl(2 * currentParameters[2]);

            try {
                out = upperTail ? Math.log(tDist.cumulativeProbability(-scaledValue))
                        : Math.log(tDist.cumulativeProbability(scaledValue));
            } catch (MathException e){
                throw new RuntimeException(e.toString());
            }

        } else {

            out =  upperTail ? NormalDistribution.standardCDF(-scaledValue, true) :
                    NormalDistribution.standardCDF(scaledValue, true);

        }
        return out;
    }


    private void update(double newData){
        dataValues.add(newData);

        double[] originalParameters=hyperprior.getParameters();
        double lambda_0 = originalParameters[1];

        double oldMu = currentParameters[0];
        double oldLambda = currentParameters[1];
        double oldAlpha = currentParameters[2];
        double oldBeta = currentParameters[3];

        double count = dataValues.size();

        double newMu = (newData - oldMu)/(lambda_0 + count) + oldMu;
        double newLambda = oldLambda + 1;
        double newAlpha = oldAlpha + 0.5;
        double newBeta = oldBeta + oldLambda*Math.pow(newData - oldMu, 2)/(2*(oldLambda+1));

        posteriorMean.setParameterValue(0, newMu);
        posteriorBeta.setParameterValue(0, newBeta);
        posteriorExpectedPrecision.setParameterValue(0, newAlpha/newBeta);

        currentParameters = new double[]{newMu, newLambda, newAlpha, newBeta};
    }



    public double calculateLogLikelihood(double[] values){

        int count = values.length;

        double[] infPredictiveDistributionParameters=hyperprior.getParameters();

        double mu_0 = infPredictiveDistributionParameters[0];
        double lambda_0 = infPredictiveDistributionParameters[1];
        double alpha_0 = infPredictiveDistributionParameters[2];
        double beta_0 = infPredictiveDistributionParameters[3];

        double lambda_n = lambda_0 + count;
        double alpha_n = alpha_0 + count/2;
        double sum = 0;
        for (Double infPeriod : values) {
            sum += infPeriod;
        }
        double mean = sum/count;

        double sumOfDifferences = 0;
        for (Double infPeriod : values) {
            sumOfDifferences += Math.pow(infPeriod-mean,2);
        }

        posteriorMean.setParameterValue(0, (lambda_0*mu_0 + sum)/(lambda_0 + count));

        double beta_n = beta_0 + 0.5*sumOfDifferences
                + lambda_0*count*Math.pow(mean-mu_0, 2)/(2*(lambda_0+count));

        posteriorBeta.setParameterValue(0, beta_n);
        posteriorExpectedPrecision.setParameterValue(0, alpha_n/beta_n);

        logL = GammaFunction.logGamma(alpha_n)
                - GammaFunction.logGamma(alpha_0)
                + alpha_0*Math.log(beta_0)
                - alpha_n*Math.log(beta_n)
                + 0.5*Math.log(lambda_0)
                - 0.5*Math.log(lambda_n)
                - (count/2)*Math.log(2*Math.PI);


        return logL;
    }

    public LogColumn[] getColumns() {
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorMean"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorMean.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorBeta"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorBeta.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorExpectedPrecision"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorExpectedPrecision.getParameterValue(0));
            }
        });

        return columns.toArray(new LogColumn[columns.size()]);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NORMAL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            double mu = xo.getDoubleAttribute(MU);
            double lambda = xo.getDoubleAttribute(LAMBDA);
            double alpha = xo.getDoubleAttribute(ALPHA);
            double beta = xo.getDoubleAttribute(BETA);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new NormalPeriodPriorDistribution(id, log, mu, lambda, alpha, beta);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false),
                AttributeRule.newDoubleRule(MU, false),
                AttributeRule.newDoubleRule(LAMBDA, false),
                AttributeRule.newDoubleRule(ALPHA, false),
                AttributeRule.newDoubleRule(BETA, false)
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles being drawn from the prior posterior distribution" +
                    "of a normal distribution of unknown mean and variance";
        }

        public Class getReturnType() {
            return NormalPeriodPriorDistribution.class;
        }
    };

}
