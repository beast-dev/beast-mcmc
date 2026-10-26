/*
 * TreeDataLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood;

import dr.evomodel.treedatalikelihood.discrete.LaplaceEvidence;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.*;
import dr.xml.*;

/**
 * The Laplace approximation of the log evidence of the maximizer's objective, re-maximized whenever the model
 * changes. The approximation is LaplaceEvidence's (plan 8 X3): at the mode of the transformed density, with the
 * Cholesky log determinant of its negative Hessian in the transformed coordinates; it replaces the earlier diagonal
 * determinant, whose sign, like that of the log-Jacobian, was wrong.
 *
 * @author Alexander Fisher
 */

public class ApproximateTreeDataLikelihood extends AbstractModelLikelihood {
    private double marginalLikelihood;
    private MaximizerWrtParameter maximizer;
    private Parameter parameter;
    private boolean likelihoodKnown = false;


    // begin parser stuff
    public static final String APPROXIMATE_LIKELIHOOD = "approximateTreeDataLikelihood";

    // end parser stuff
    public ApproximateTreeDataLikelihood(MaximizerWrtParameter maximizer) {
        super(APPROXIMATE_LIKELIHOOD);

        if (maximizer.getTransform() != null && !maximizer.includesJacobian()) {
            throw new IllegalArgumentException("The Laplace approximation needs the mode of the transformed density: " +
                    "give the maximizer includeJacobian=\"true\"");
        }
        this.maximizer = maximizer;
        this.parameter = maximizer.getGradient().getParameter();
        updateParameterMAP();
        updateMarginalLikelihood();
        addVariable(parameter);
        if (maximizer.getLikelihood().getModel() != null) {
            // a density of the parameter alone (a distributionLikelihood) has no model
            addModel(maximizer.getLikelihood().getModel());
        }
    }

    private void updateMarginalLikelihood() {
        this.marginalLikelihood = new LaplaceEvidence(maximizer).getLogEvidence();
        likelihoodKnown = true;
    }

    private void updateParameterMAP() {
        maximizer.maximize();
    }


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return APPROXIMATE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MaximizerWrtParameter maximizer =
                    (MaximizerWrtParameter) xo.getChild(MaximizerWrtParameter.class);

            try {
                return new ApproximateTreeDataLikelihood(maximizer);
            } catch (IllegalArgumentException | IllegalStateException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Approximates the marginal likelihood of the data given the tree using Laplace approximation";
        }

        @SuppressWarnings("rawtypes")
        @Override
        public Class getReturnType() {
            return ApproximateTreeDataLikelihood.class;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MaximizerWrtParameter.class)
        };
    };

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            updateParameterMAP();
            updateMarginalLikelihood();
        }
        return marginalLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }
}