/*
 * IntegratedFactorTraitDataModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class IntegratedFactorAnalysisLikelihood extends AbstractModelLikelihood implements ContinuousTraitPartialsProvider {

    public IntegratedFactorAnalysisLikelihood(String name,
                                              MatrixParameterInterface loadings,
                                              Parameter traitPrecision) {
        super(name);

        this.loadings = loadings;
        this.traitPrecision = traitPrecision;
        this.missingIndices = null; // TODO

        this.numTraits = 0; // TODO
        this.dimTrait = 0; // TODO
    }

    @Override
    public boolean bufferTips() {
        return true;
    }

    @Override
    public int getTraitCount() {
        return numTraits;
    }

    @Override
    public int getTraitDimension() {
        return dimTrait;
    }

    @Override
    public PrecisionType getPrecisionType() {
        return PrecisionType.FULL;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {
        return new double[0];
    }

    @Override
    public double[] getTipPartial(int taxonIndex) {
        return new double[0];
    }

    @Override
    public double[] getTipObservation(int taxonIndex, PrecisionType precisionType) {
        return new double[0];
    }

    @Override
    public List<Integer> getMissingIndices() {
        return missingIndices; // TODO Fix use-case (all values should be missing)
    }

    @Override
    public CompoundParameter getParameter() {
        return null; // TODO Fix use-case
    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood() {
        return 0.0;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == loadings || variable == traitPrecision) {
            fireModelChanged(this);
//            fireModelChanged(this, getTaxonIndex(index));
        } else {
            throw new RuntimeException("Unhandled parameter change type");
        }
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnow = likelihoodKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnow;
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnow;

    private final int numTraits;
    private final int dimTrait;
    private final MatrixParameterInterface loadings;
    private final Parameter traitPrecision;
    private final List<Integer> missingIndices;

    // TODO Move remainder into separate class file
    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getElementFirstChild(LOADINGS);
            Parameter traitPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

            return new IntegratedFactorAnalysisLikelihood(xo.getId(), loadings, traitPrecision);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return IntegratedFactorAnalysisLikelihood.class;
        }

        @Override
        public String getParserName() {
            return INTEGRATED_FACTOR_Model;
        }
    };

    public static final String INTEGRATED_FACTOR_Model = "integratedFactorModel";
    public static final String LOADINGS = "loadings";
    public static final String PRECISION = "precision";

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(LOADINGS, new XMLSyntaxRule[] {
                    new ElementRule(MatrixParameterInterface.class),
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),


    };
}
