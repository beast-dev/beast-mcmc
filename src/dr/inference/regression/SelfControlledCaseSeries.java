/*
 * SelfControlledCaseSeries.java
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

package dr.inference.regression;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.*;

/**
 * @author Marc Suchard
 * @author Trevor Shaddox
 */
public class SelfControlledCaseSeries extends AbstractModelLikelihood {

    public static final String SCCS_NAME = "selfControlledCaseSeries";
    public static final String FILE_NAME = "fileName";
    public static final String BETA = "beta";
    public static final String PRECISION = "precision";

    public SelfControlledCaseSeries(String name, String fileName, Parameter beta, Parameter precision) {
        super(name);

        regressionInterface = RegressionJNIWrapper.loadLibrary();

        // Load data and find mode
        instance = regressionInterface.loadData(fileName);
        regressionInterface.setPriorType(instance, RegressionJNIWrapper.NORMAL_PRIOR);

        this.precision = precision;
        setPrecision();
        precisionChanged = true;

        // Set beta to mode
        final int dim = regressionInterface.getBetaSize(instance);
        if (dim != beta.getDimension()) {
            beta.setDimension(dim);
        }

        // Start beta at mode (given precision)
        this.beta = beta;
        double[] mode = getMode();
        for (int i = 0; i < beta.getDimension(); ++i) {
            beta.setParameterValue(i, mode[i]);
        }
        logSCCSLikelihood = regressionInterface.getLogLikelihood(instance);
        logSCCSPrior = regressionInterface.getLogPrior(instance);
        betaChanged = false; // Internal state is at mode

        addVariable(beta);
        addVariable(precision);
    }

    private void setPrecision() {
        regressionInterface.setHyperprior(instance, 1.0 / precision.getParameterValue(0));
    }

    public double[] getMode() {
        if (precisionChanged) {
            setPrecision();
            mode = null;
        }
        if (mode == null) {
            regressionInterface.findMode(instance);
            mode = new double[beta.getDimension()];
            for (int i = 0; i < beta.getDimension(); ++i) {
                mode[i] = regressionInterface.getBeta(instance, i);
            }
            betaChanged = true; // Internal beta-state is at mode, not betaParameter
//            betaFlag.clear();
            newMode = true;
//            System.err.println("A");
            if (DEBUG_MODE) {
                System.err.println("Recomputed mode!");
            }
        }
        double[] rtn = new double[mode.length];
        System.arraycopy(mode, 0, rtn, 0, mode.length);
        return rtn;
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
        if (variable == beta) {
            betaChanged = true;
            if (type == Variable.ChangeType.ALL_VALUES_CHANGED) {
                betaFlag.clear();
            } else {
                betaFlag.add(index);
            }
        } else if (variable == precision) {
            precisionChanged = true;
        } else {
            throw new IllegalArgumentException("Unknown variable in SCCS");
        }
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        storedLogSCCSLikelihood = logSCCSLikelihood;
        storedLogSCCSPrior = logSCCSPrior;

        storedBetaChanged = betaChanged;
        storedPrecisionChanged = precisionChanged;

        storedPrecision = precision.getParameterValue(0);
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        logSCCSLikelihood = storedLogSCCSLikelihood;
        logSCCSPrior = storedLogSCCSPrior;

        betaChanged = storedBetaChanged;
        precisionChanged = storedPrecisionChanged;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {
        if (storedPrecision != precision.getParameterValue(0)) {
            mode = null; // Accepted new precision state; mode has moved
        }
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
        return calculateLogLikelihood();
    }

    private double calculateLogLikelihood() {

        if (betaChanged) {
            if (betaFlag.isEmpty() || newMode) {
                regressionInterface.setBeta(instance, beta.getParameterValues());
                newMode = false;
            } else {
                while (!betaFlag.isEmpty()) {
                    final int index = betaFlag.remove();
                    regressionInterface.setBeta(instance, index, beta.getParameterValue(index));
                }
            }
        }

        if (precisionChanged) {
            setPrecision();
        }

        if (betaChanged) {
            logSCCSLikelihood = regressionInterface.getLogLikelihood(instance);
        }

        if (betaChanged || precisionChanged) {
            logSCCSPrior = regressionInterface.getLogPrior(instance);
        }

        betaChanged = false;
        precisionChanged = false;

        double logLike = logSCCSLikelihood + logSCCSPrior;

        if (DEBUG_LAZY) {
            double checkLike = regressionInterface.getLogLikelihood(instance);
            double checkPrior = regressionInterface.getLogPrior(instance);

            double check = checkLike + checkPrior;
            if (check != logLike) {
                System.err.println("Error in internal state in calculateLogLikelihood()");
                System.err.println(checkLike + " " + logSCCSLikelihood + " d: " + (checkLike - logSCCSLikelihood));
                System.err.println(checkPrior + " " + logSCCSPrior + " d: " + (checkPrior - logSCCSPrior));
                System.err.println(betaChanged + " " + precisionChanged);
                System.exit(-1);
            }
        }

        return logLike;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        betaChanged = true;
        newMode = true;
        precisionChanged = true;
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

    private double logSCCSLikelihood;
    private double logSCCSPrior;
    private double storedLogSCCSLikelihood;
    private double storedLogSCCSPrior;

    private boolean betaChanged;
    private boolean precisionChanged;
    private boolean storedBetaChanged;
    private boolean storedPrecisionChanged;

//    private boolean[] betaFlag;
//    private boolean betaFlagAll;
    private Queue<Integer> betaFlag = new LinkedList<Integer>();
    private boolean newMode = false;

    private double[] mode = null;
    private double storedPrecision;

    private static final boolean DEBUG_MODE = false;
    private static final boolean DEBUG_LAZY = false;

}
