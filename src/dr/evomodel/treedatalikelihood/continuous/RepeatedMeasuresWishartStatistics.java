/*
 * RepeatedMeasuresWishartStatistics.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresWishartStatistics implements ConjugateWishartStatisticsProvider {


    private final FullPrecisionContinuousTraitPartialsProvider traitModel;
    private final ConditionalTraitSimulationHelper extensionHelper;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final double[] outerProduct;
    private final int dimTrait;
    private double[] buffer;
    private boolean forceResample;

    public RepeatedMeasuresWishartStatistics(FullPrecisionContinuousTraitPartialsProvider traitModel,
                                             TreeDataLikelihood treeLikelihood,
                                             boolean forceResample) {
        this.traitModel = traitModel;

        this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();
        this.extensionHelper = likelihoodDelegate.getExtensionHelper();

        this.dimTrait = traitModel.getTraitDimension();

        this.outerProduct = new double[dimTrait * dimTrait];

        this.forceResample = forceResample;

    }


    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return traitModel.getExtensionPrecisionParameter();
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {

        if (forceResample) {
            likelihoodDelegate.fireModelChanged();
        }

        ConditionalTraitSimulationHelper.JointSamples traits = extensionHelper.drawTraitsAboveAndBelow(traitModel, true);

        double[] valuesAbove = traits.getTraitsAbove();
        double[] valuesBelow = traits.getTraitsBelow();

        int nTipsTotal = valuesAbove.length / dimTrait;

        if (buffer == null) {
            buffer = new double[dimTrait * nTipsTotal];
        }


        DenseMatrix64F XminusY = DenseMatrix64F.wrap(nTipsTotal, dimTrait, buffer);
        DenseMatrix64F X = DenseMatrix64F.wrap(nTipsTotal, dimTrait, valuesAbove);
        DenseMatrix64F Y = DenseMatrix64F.wrap(nTipsTotal, dimTrait, valuesBelow);

        CommonOps.subtract(X, Y, XminusY);

        DenseMatrix64F outerProductMat = DenseMatrix64F.wrap(dimTrait, dimTrait, outerProduct);

        CommonOps.multTransA(XminusY, XminusY, outerProductMat);


        return new WishartSufficientStatistics(nTipsTotal, outerProduct);
    }

    public void setForceResample(Boolean b) {
        forceResample = b;
    }


    private static final boolean DEBUG = false;


    public static final String RM_WISHART_STATISTICS = "repeatedMeasuresWishartStatistics";
    private static final String FORCE_RESAMPLE = "forceResample";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            RepeatedMeasuresTraitDataModel traitModel =
                    (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);

            boolean forceResample = xo.getAttribute(FORCE_RESAMPLE, true);

            return new RepeatedMeasuresWishartStatistics(traitModel, dataLikelihood, forceResample);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newBooleanRule(FORCE_RESAMPLE, true)
        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return null;
        }

        @Override
        public String getParserName() {
            return RM_WISHART_STATISTICS;
        }
    };
}


