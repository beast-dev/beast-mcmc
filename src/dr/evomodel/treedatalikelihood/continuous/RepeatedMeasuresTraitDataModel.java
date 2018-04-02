/*
 * ContinuousTraitData.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresTraitDataModel extends
        ContinuousTraitDataModel implements ContinuousTraitPartialsProvider {

    private final String traitName;
    private final Parameter samplingPrecision;

    public RepeatedMeasuresTraitDataModel(String name,
                                          CompoundParameter parameter,
                                          List<Integer> missingIndices,
                                          boolean useMissingIndices,
                                          final int dimTrait,
                                          Parameter samplingPrecision) {
        super(name, parameter, missingIndices, useMissingIndices, dimTrait, PrecisionType.FULL);
        this.traitName = name;
        this.samplingPrecision = samplingPrecision;
        addVariable(samplingPrecision);

        samplingPrecision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                samplingPrecision.getDimension()));

        if (samplingPrecision.getDimension() != dimTrait) {
            throw new RuntimeException("Currently only implemented for diagonal deflation");
        }
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        assert (numTraits == 1);

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);
        if (fullyObserved == true){
            return partial;
        }
        DenseMatrix64F V = MissingOps.wrap(partial,dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

        for (int index = 0; index< dimTrait; index++){
            V.set(index, index, V.get(index, index) + 1 / samplingPrecision.getParameterValue(index));
        }

        DenseMatrix64F P = new DenseMatrix64F(dimTrait, dimTrait);
        MissingOps.safeInvert(V, P, false);

        MissingOps.unwrap(P, partial, dimTrait);
        MissingOps.unwrap(V, partial, dimTrait + dimTrait * dimTrait);

        return partial;
    }

    public Parameter getSamplingPrecision() { return samplingPrecision; }

    public String getTraitName() { return traitName; }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);

        if (variable == samplingPrecision) {
            fireModelChanged();
        }
    }
    
    // TODO Move remainder into separate class file
    private static final String REPEATED_MEASURES_MODEL = "repeatedMeasuresModel";
    private static final String PRECISION = "samplingPrecision";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(TreeModel.class);
            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME,
                            treeModel, true);
            CompoundParameter traitParameter = returnValue.traitParameter;
            List<Integer> missingIndices = returnValue.missingIndices;

            Parameter samplingPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

            String traitName = returnValue.traitName;
            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
                    xo.getChild(MultivariateDiffusionModel.class);

            return new RepeatedMeasuresTraitDataModel(
                    traitName,
                    traitParameter,
                    missingIndices,
                    true,
                    diffusionModel.getPrecisionParameter().getRowDimension(),
                    samplingPrecision
            );
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
            return RepeatedMeasuresTraitLikelihood.class;
        }

        @Override
        public String getParserName() {
            return REPEATED_MEASURES_MODEL;
        }
    };

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(PRECISION, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            // Tree trait parser
            new ElementRule(MutableTreeModel.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(MultivariateDiffusionModel.class),
    };
}
