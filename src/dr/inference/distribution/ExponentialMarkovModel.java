/*
 * ExponentialMarkovModel.java
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

package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that acts as a model for exponentially distributed data.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Gerton Lunter
 * @version $Id: ExponentialMarkovModel.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */

public class ExponentialMarkovModel extends AbstractModel implements Likelihood {

    public static final String EXPONENTIAL_MARKOV_MODEL = "exponentialMarkovLikelihood";
    public static final String CHAIN_PARAMETER = "chainParameter";
    public static final String JEFFREYS = "jeffreys";
    public static final String REVERSE = "reverse";
    public static final String SHAPE = "shape";

    /**
     * Constructor.
     */
    public ExponentialMarkovModel(Parameter chainParameter, boolean jeffreys, boolean reverse, double shape) {

        super(EXPONENTIAL_MARKOV_MODEL);

        this.chainParameter = chainParameter;
        this.jeffreys = jeffreys;
        this.reverse = reverse;
        this.shape = shape;

        addParameter(chainParameter);
        chainParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, chainParameter.getDimension()));
    }

    public Parameter getChainParameter() {
        return getParameter(0);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EXPONENTIAL_MARKOV_MODEL;
        }

        /**
         * Reads a gamma distribution model from a DOM Document element.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject object = (XMLObject) xo.getChild(CHAIN_PARAMETER);

            Parameter chainParameter = (Parameter) object.getChild(0);

            boolean jeffreys = xo.getAttribute(JEFFREYS, false);
            boolean reverse = xo.getAttribute(REVERSE, false);

            double shape = xo.getAttribute(SHAPE, 1.0);
            if (shape < 1.0) {
                throw new XMLParseException("ExponentialMarkovModel: shape parameter must be >= 1.0");
            }

            if (shape == 1.0) {
                System.out.println("Exponential markov model on parameter " + chainParameter.getParameterName() + " (jeffreys=" + jeffreys + ", reverse=" + reverse + ")");
            } else {
                System.out.println("Gamma markov model on parameter " + chainParameter.getParameterName() + " (jeffreys=" + jeffreys + ", reverse=" + reverse + " shape=" + shape + ")");
            }

            return new ExponentialMarkovModel(chainParameter, jeffreys, reverse, shape);
        }

        public String getParserDescription() {
            return "A continuous state, discrete time markov chain in which each new state is an exponentially distributed variable with a mean of the previous state.";
        }

        public Class getReturnType() {
            return ExponentialMarkovModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(JEFFREYS, true),
                AttributeRule.newBooleanRule(REVERSE, true),
                new ElementRule(CHAIN_PARAMETER, Parameter.class)
        };
    };

    // **************************************************************
    // Likelihood
    // **************************************************************

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;
    }

    private int index(int i) {
        if (reverse)
            return chainParameter.getDimension() - i - 1;
        else
            return i;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {

        double logL = 0.0;
        // jeffreys Prior!
        if (jeffreys) {
            logL += -Math.log(chainParameter.getParameterValue(index(0)));
        }
        for (int i = 1; i < chainParameter.getDimension(); i++) {
            final double mean = chainParameter.getParameterValue(index(i - 1));
            final double x = chainParameter.getParameterValue(index(i));
            //logL += dr.math.distributions.ExponentialDistribution.logPdf(x, 1.0/mean);

            final double scale = mean / shape;
            logL += GammaDistribution.logPdf(x, shape, scale);
        }
        return logL;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private Parameter chainParameter = null;
    private boolean jeffreys = false;
    private boolean reverse = false;
    private double shape = 1.0;

}

