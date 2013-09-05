/*
 * SelfControlledCaseSeries.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.regression;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Trevor Shaddox
 */
public class SelfControlledCaseSeries extends AbstractModelLikelihood {

    public static final String SCCS_NAME = "selfControlledCaseSeries";
    public static final String FILE_NAME = "fileName";
    public static final String BETA = "beta";
    public static final String PRECISION = "precision";

    /**
     * @param name Model Name
     * @param fileName
     * @param beta
     * @param precision
     */
    public SelfControlledCaseSeries(String name, String fileName, Parameter beta, Parameter precision) {
        super(name);

        regressionInterface = RegressionJNIWrapper.loadLibrary();

        // Load data and find mode
        instance = regressionInterface.loadData(fileName);
        regressionInterface.setPriorType(instance, RegressionJNIWrapper.NORMAL_PRIOR);
        regressionInterface.setHyperprior(instance, 1.0 / precision.getParameterValue(0));
        this.precision = precision;
        precisionChanged = false;

        regressionInterface.findMode(instance);

        // Set beta to mode
        final int dim = regressionInterface.getBetaSize(instance);
        if (dim != beta.getDimension()) {
            beta.setDimension(dim);
        }
        for (int i = 0; i < beta.getDimension(); ++i) {
            beta.setParameterValue(i, regressionInterface.getBeta(instance, i));
        }
        this.beta = beta;
        betaChanged = false;

        addVariable(beta);
        addVariable(precision);
        logLikelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        logLikelihoodKnown = false;
        betaChanged = true;
        precisionChanged = true;
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        savedLogLikelihoodKnown = logLikelihoodKnown;
        savedLogLikelihood = logLikelihood;
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        logLikelihoodKnown = false;

//        logLikelihoodKnown = savedLogLikelihoodKnown;
//        logLikelihood = savedLogLikelihood;
        betaChanged = true;
        precisionChanged = true;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {
        // Do nothing
    }

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {
        if (!logLikelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            logLikelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {
        // TODO Check is beta has changed
        for (int i = 0; i < beta.getDimension(); ++i) {
            regressionInterface.setBeta(instance, i, beta.getParameterValue(i));
        }
        // TODO Check if precision has changed
        regressionInterface.setHyperprior(instance, 1.0 / precision.getParameterValue(0));

        return regressionInterface.getLogLikelihood(instance) + regressionInterface.getLogPrior(instance);
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        logLikelihoodKnown = false;
        regressionInterface.makeDirty(instance);
    }


  public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SCCS_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            Parameter beta = (Parameter) xo.getElementFirstChild(BETA);
            Parameter precision = (Parameter) xo.getElementFirstChild(PRECISION);

            return new SelfControlledCaseSeries(xo.getId(), fileName, beta, precision);
        }

        public String getParserDescription() {
            return "Self-controlled case series design.";
        }

        public Class getReturnType() {
            return SelfControlledCaseSeries.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME),
                new ElementRule(BETA, Parameter.class),
                new ElementRule(PRECISION, Parameter.class),
        };
    };

    private final RegressionJNIWrapper regressionInterface;
    private final int instance;
    private final Parameter beta;
    private final Parameter precision;

    private boolean logLikelihoodKnown;
    private boolean savedLogLikelihoodKnown;
    private double logLikelihood;
    private double savedLogLikelihood;

    private boolean betaChanged;
    private boolean precisionChanged;
}
