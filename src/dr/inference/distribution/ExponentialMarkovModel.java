/*
 * ExponentialMarkovModel.java
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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GammaDistribution;
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

public class ExponentialMarkovModel extends AbstractModelLikelihood {

    public static final String EXPONENTIAL_MARKOV_MODEL = "exponentialMarkovLikelihood";

    /**
     * Constructor.
     */
    public ExponentialMarkovModel(Parameter chainParameter, boolean jeffreys, boolean reverse, double shape) {

        super(EXPONENTIAL_MARKOV_MODEL);

        this.chainParameter = chainParameter;
        this.jeffreys = jeffreys;
        this.reverse = reverse;
        this.shape = shape;

        addVariable(chainParameter);
        chainParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, chainParameter.getDimension()));
    }

    public Parameter getChainParameter() {
        return (Parameter)getVariable(0);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
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

