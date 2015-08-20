/*
 * MultivariatePolyaDistributionLikelihood.java
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

package dr.math.distributions;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.GammaFunction;
import dr.xml.*;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Package: MultivariatePolyaDistributionLikelihood
 * Description:
 * this class provides a model for over-dispersed multinomial counts. The model follows Dirichlet-Multinomial distribution with
 * multinomial parameters integrated out analytically. This model is also known as Multivariate Polya distribution.
 * Standard parametrization involves k intensities a_i's. This implementation uses the standard parametrization internally, but allows
 * for re-parametrization as frequencies (k-1 df) and dispersion parameters, where a = sum_i=1^k a_i is dispersion and f_i = a_i/a
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Sep 23, 2010
 *         Time: 12:25:14 PM
 */
public class MultivariatePolyaDistributionLikelihood extends AbstractModel implements Likelihood {

    protected Parameter frequencies;
    protected Parameter dispersion;
    protected Parameter alphas;
    protected boolean usingAlphas;
    protected boolean isAlphasKnown;
    protected MatrixParameter data;
    protected double fixedNorm;
    protected double variableNorm;
    protected double storedFixedNorm;
    protected double storedVariableNorm;
    protected double logLikelihood;
    protected double storedLogLikelihood;
    protected boolean isLogLikelihoodKnown;
    protected boolean isFixedNormKnown;
    protected boolean isVariableNormKnown;
    protected double rowSums[];


    public MultivariatePolyaDistributionLikelihood(String modelID, MatrixParameter data, Parameter frequencies, Parameter dispersion) {
        super(modelID);
        this.frequencies = frequencies;
        this.dispersion = dispersion;
        this.alphas = new Parameter.Default(frequencies.getDimension());
        computeAlphas();
        this.data = data;
        isFixedNormKnown = false;
        isVariableNormKnown = false;
        addVariable(this.frequencies);
        addVariable(this.dispersion);
        addVariable(this.data);
        if (this.alphas.getDimension() != data.getColumnDimension()) {
            System.err.println("Dimensions of the frequency vector and number of columns do not match!");
        }
    }

    public MultivariatePolyaDistributionLikelihood(String modelID, MatrixParameter data, Parameter alphas) {
        super(modelID);
        this.alphas = alphas;
        isAlphasKnown = true;
        usingAlphas = true;

        this.frequencies = new Parameter.Default(alphas.getDimension());
        this.dispersion = new Parameter.Default(1);
        this.data = data;
        isFixedNormKnown = false;
        isVariableNormKnown = false;
        addVariable(this.alphas);
        addVariable(this.data);
        if (this.alphas.getDimension() != data.getColumnDimension()) {
            System.err.println("Dimensions of the frequency vector and number of columns do not match!");
        }
    }

    /* Compute alphas from frequencies and dispersion
    */
    protected void computeAlphas(){
        double disp=dispersion.getParameterValue(0);
        double[] freqs = frequencies.getParameterValues();

        for(int i=0; i<alphas.getDimension(); ++i){
            alphas.setParameterValueQuietly(i, disp*freqs[i]);
        }
        alphas.setParameterValueNotifyChangedAll(0, alphas.getParameterValue(0));
        isAlphasKnown = true;
    }

    public MultivariatePolyaDistributionLikelihood(String modelID) {
        super(modelID);
    }

    public double calculateLogLikelihood() {
        // R code for this function:
        // //assuming X[,1] is row totals, lfactX is log factorial of X //
        // logLikes = sapply(1:n, function(i)  lfactX[subset,1][i] - sum(lfactX[subset,][i,2:p1]) +
        // lgamma(sum(alpha)) - lgamma(X[subset,1][i] + sum(alpha)) + sum(lgamma(X[subset,2:p1][i,] + alpha)) - sum(lgamma(alpha)));
        // sum(logLikes)
        if(!isAlphasKnown) computeAlphas();
        if (!isFixedNormKnown) {
            computeFixedNorm();
        }
        if (!isVariableNormKnown) {
            computeVariableNorm();
        }
        double logP = fixedNorm + variableNorm;

        double disp = 0;
        double[] a = alphas.getParameterValues();
        for(int i = 0; i< alphas.getDimension(); ++i){
            disp = disp + a[i];
        }

        
        for (int i = 0; i < data.getRowDimension(); ++i) {
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                logP += GammaFunction.lnGamma(data.getParameterValue(i, j) + a[j]);
            }
            logP -= GammaFunction.lnGamma(rowSums[i] + disp);
        }
        return logP;
    }

    protected void computeFixedNorm() {
        rowSums = new double[data.getRowDimension()];
        for (int i = 0; i < data.getRowDimension(); ++i) {
            rowSums[i] = 0;
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                rowSums[i] += data.getParameterValue(i, j);
            }
        }

        fixedNorm = 0;
        for (int i = 0; i < data.getRowDimension(); ++i) {
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                fixedNorm -= GammaFunction.lnGamma(data.getParameterValue(i, j) + 1);
            }
            fixedNorm += GammaFunction.lnGamma(rowSums[i] + 1);
        }

        isFixedNormKnown = true;
    }

    protected void computeVariableNorm() {
        double disp = 0;
        double[] a = alphas.getParameterValues();
        for(int i = 0; i< alphas.getDimension(); ++i){
            disp = disp + a[i];
        }
        variableNorm = GammaFunction.lnGamma(disp);

        for (int i = 0; i < alphas.getDimension(); ++i) {
            variableNorm -= GammaFunction.lnGamma(a[i]);
        }
        variableNorm *= data.getRowDimension();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable.getVariableName().equals(frequencies.getVariableName()) || variable.getVariableName().equals(dispersion.getVariableName())) {
            isAlphasKnown = false;
            isVariableNormKnown = false;
        } else if (variable.getVariableName().equals(data.getVariableName())) {
            isFixedNormKnown = false;
        }
        else if(variable.getVariableName().equals(alphas.getVariableName())){
            isVariableNormKnown = false;
        }
    }

    protected void storeState() {
        storedVariableNorm = variableNorm;
        storedFixedNorm = fixedNorm;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        variableNorm = storedVariableNorm;
        fixedNorm = storedFixedNorm;
        logLikelihood = storedLogLikelihood;
        if(!usingAlphas) computeAlphas();
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!isLogLikelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
        }
        return logLikelihood;
    }

    public void makeDirty() {
        isLogLikelihoodKnown = false;
        isVariableNormKnown = false;
        isFixedNormKnown = false;
    }

    public String prettyName() {
        return "Multivariate Polya Distribution Likelihood";
    }

    public boolean evaluateEarly() {
        return false;
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
    }

    public void setUsed() {
    }


    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new NumberColumn(this.getId()) {
                    public double getDoubleValue() {
                        return getLogLikelihood();
                    }
                }
        };
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MVPLIKE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameter data;
            Parameter dispersion;
            Parameter frequencies;
            Parameter rates;

            if (xo.hasChildNamed(DATA)) {
                data = (MatrixParameter) xo.getChild(DATA).getChild(MatrixParameter.class);
            } else {
                throw new XMLParseException("Missing data element!");
            }

            if (xo.hasChildNamed(RATES)) {
                rates = (Parameter) xo.getChild(RATES).getChild(Parameter.class);
                if (rates.getDimension() != data.getColumnDimension()) {
                    throw new XMLParseException("The number of data columns must match the dimension of " + RATES
                            + " parameter (" + data.getColumnDimension() + " != " + rates.getDimension() + "!");
                }
            }
            else if (xo.hasChildNamed(FREQ)) {
                frequencies = (Parameter) xo.getChild(FREQ).getChild(Parameter.class);
                if (xo.hasChildNamed(DISPERSION)) {
                    dispersion = (Parameter) xo.getChild(DISPERSION).getChild(Parameter.class);
                } else {
                    throw new XMLParseException(DISPERSION + " element has to be specified when using " + FREQ
                            +" parametrization");
                }
                if (dispersion.getDimension() != 1) {
                    throw new XMLParseException("Dispersion parameter must be of dimmension exactly 1!");
                }

                if (frequencies.getDimension() != data.getColumnDimension()) {
                    throw new XMLParseException("The number of data columns must match the dimension of "+ FREQ
                            + " parameter (" + data.getColumnDimension() + " != " + frequencies.getDimension() + "!");
                }
                return new MultivariatePolyaDistributionLikelihood(MVPLIKE, data, frequencies, dispersion);

            } else {
                throw new XMLParseException("Either " + FREQ + " or " + RATES + "element has to be specified!");
            }


            return new MultivariatePolyaDistributionLikelihood(MVPLIKE, data, rates);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, false),
                new XORRule(new ElementRule(RATES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, false),
                        new ElementRule(FREQ, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, false)),
                new ElementRule(DISPERSION, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
        };

        public Class getReturnType() {
            return MatrixParameter.class;
        }
    };
    public static final String MVPLIKE = "mvPolyaLikelihood";
    public static final String DATA = "data";
    public static final String DISPERSION = "dispersion";
    public static final String FREQ = "frequencies";
    public static final String RATES = "alpha";
}
