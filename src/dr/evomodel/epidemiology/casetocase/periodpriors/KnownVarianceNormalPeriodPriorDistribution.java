/*
 * KnownVarianceNormalPeriodPriorDistribution.java
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
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

/**
 The assumption here is that the periods are drawn from a normal distribution with unknown mean and variance.
 The hyperprior is the conjugate, normal-gamma distribution.

 @author Matthew Hall
 */
public class KnownVarianceNormalPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String NORMAL = "knownVarianceNormalPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    // This gets confusing. The data is assumed to be normally distributed with mean mu and stdev sigma. Sigma is
    // known. The prior on mu is that it is _also_ normally distributed with mean mu_0 and stdev sigma_0.

    public static final String MU_0 = "mu0";
    public static final String SIGMA = "sigma";
    public static final String SIGMA_0 = "sigma0";

    private NormalDistribution hyperprior;

    private Parameter posteriorMean;
    private Parameter posteriorVariance;
    private double sigma;

    private ArrayList<Double> dataValues;
    private double[] currentParameters;

    public KnownVarianceNormalPeriodPriorDistribution(String name, boolean log, double sigma,
                                                      NormalDistribution hyperprior){
        super(name, log);
        this.hyperprior = hyperprior;
        posteriorVariance = new Parameter.Default(1);
        posteriorMean = new Parameter.Default(1);
        addVariable(posteriorVariance);
        addVariable(posteriorMean);
        this.sigma = sigma;
    }

    public KnownVarianceNormalPeriodPriorDistribution(String name, boolean log, double sigma,
                                                      double mu_0, double sigma_0){
        this(name, log, sigma, new NormalDistribution(mu_0, sigma_0));
    }

    public void reset(){
        dataValues = new ArrayList<Double>();
        currentParameters[0] = hyperprior.getMean();
        currentParameters[1] = hyperprior.getSD();
        logL = 0;
    }

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
        double sd = currentParameters[1];

        return NormalDistribution.logPdf(value, mean, Math.sqrt(Math.pow(sd, 2) + Math.pow(sigma, 2)));
    }

    public double calculateLogPosteriorPredictiveCDF(double value, boolean upperTail){
        double mean = currentParameters[0];
        double sd = currentParameters[1];

        double scaledValue = (value - mean)/Math.sqrt(Math.pow(sd, 2) + Math.pow(sigma, 2));

        return upperTail ?  NormalDistribution.standardCDF(-scaledValue, true) :
                NormalDistribution.standardCDF(scaledValue, true);
    }

    private void update(double newData){
        dataValues.add(newData);

        double originalMean = hyperprior.getMean();
        double originalSD = hyperprior.getSD();

        double count = dataValues.size();

        double dataMean = 0;
        for(double value: dataValues){
            dataMean += value;
        }
        dataMean /= count;

        double newSD = Math.sqrt(1/(count/Math.pow(sigma,2) + 1/Math.pow(originalSD,2)));

        double newMean = Math.pow(newSD,2)*(originalMean/Math.pow(originalSD,2) + count*dataMean/Math.pow(sigma,2));

        currentParameters = new double[]{newMean, newSD};
    }


    public double calculateLogLikelihood(double[] values){

        int count = values.length;

        double mu_0 = hyperprior.getMean();
        double sigma_0 = hyperprior.getSD();

        double var = Math.pow(sigma, 2);
        double var_0 = Math.pow(sigma_0, 2);

        double sum = 0;
        double sumOfSquares = 0;
        for (Double infPeriod : values) {
            sum += infPeriod;
            sumOfSquares += Math.pow(infPeriod, 2);
        }
        double mean = sum/count;

        posteriorMean.setParameterValue(0, ((mu_0/var_0) + sum/var)/(1/var_0 + count/var));
        posteriorVariance.setParameterValue(0, 1/(1/var_0 + count/var));

        logL = Math.log(sigma)
                - count * Math.log(Math.sqrt(2*Math.PI)*sigma)
                - Math.log(Math.sqrt(count*var_0 + var))
                + -sumOfSquares/(2*var) - Math.pow(mu_0, 2)/(2*var_0)
                + (Math.pow(sigma_0*count*mean/sigma, 2) + Math.pow(sigma*mu_0/sigma_0, 2) + 2*count*mean*mu_0)
                /(2*(count*var_0 + var));


        return logL;
    }

    public LogColumn[] getColumns() {
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorMean"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorMean.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorVariance"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorVariance.getParameterValue(0));
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

            double mu_0 = xo.getDoubleAttribute(MU_0);
            double sigma = xo.getDoubleAttribute(SIGMA);
            double sigma_0 = xo.getDoubleAttribute(SIGMA_0);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new KnownVarianceNormalPeriodPriorDistribution(id, log, sigma, mu_0, sigma_0);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false),
                AttributeRule.newDoubleRule(MU_0, false),
                AttributeRule.newDoubleRule(SIGMA, false),
                AttributeRule.newDoubleRule(SIGMA_0, false),
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles being drawn from the prior posterior distribution" +
                    "of a normal distribution of unknown mean and known standard deviation sigma";
        }

        public Class getReturnType() {
            return KnownVarianceNormalPeriodPriorDistribution.class;
        }
    };

}
