/*
 * WishartStatisticsWrapper.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.*;

import java.util.List;

import static dr.evomodel.treedatalikelihood.ProcessSimulationDelegate.ConditionalOnPartiallyMissingTipsDelegate.getPartiallyMissingTraitName;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 */
public class WishartStatisticsWrapper extends AbstractModel implements ConjugateWishartStatisticsProvider {

    public static final String PARSER_NAME = "wishartStatistics";
    public static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    private final LikelihoodTreeTraversal treeTraversalDelegate;
    private final ContinuousRateTransformation rateTransformation;

    private final int dimTrait;
    private final int numTrait;
    private final int dimPartial;
    private final int nodeCount;

    public WishartStatisticsWrapper(final String name,
                                    final String traitName,
                                    final TreeDataLikelihood dataLikelihood,
                                    final ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(name);
        this.dataLikelihood = dataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;
        this.rateTransformation = likelihoodDelegate.getRateTransformation();

        this.dimTrait = likelihoodDelegate.getTraitDim();
        this.numTrait = likelihoodDelegate.getTraitCount();
        this.nodeCount = dataLikelihood.getTree().getNodeCount();
        this.dimPartial = dimTrait * numTrait;

        addModel(dataLikelihood);

        String partialTraitName = getPartiallyMissingTraitName(traitName);
        TreeTrait trait = dataLikelihood.getTreeTrait(partialTraitName);

        hasPartiallyMissingTraits = trait != null;
        System.err.println("ps ? " + (trait instanceof ProcessSimulation ? "yes" : "no"));
        System.err.println("partial trait is " + (trait == null ? "null" : "notnull"));

        treeTraversalDelegate = new LikelihoodTreeTraversal(
                dataLikelihood.getTree(),
                dataLikelihood.getBranchRateModel(),
                TreeTraversal.TraversalType.POST_ORDER);

        if (likelihoodDelegate.getIntegrator() instanceof ContinuousDiffusionIntegrator.Multivariate) {
            outerProductDelegate = likelihoodDelegate.createObservedDataOnly(likelihoodDelegate);
        } else {
            outerProductDelegate = likelihoodDelegate;
        }

        traitDataKnown = false;
        outerProductsKnown = false;
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {

        if (!outerProductsKnown) {
            computeOuterProducts();
            outerProductsKnown = true;
        }

        return wishartStatistics;
    }

    private void simulateMissingTraits() {
        System.err.println("Attempting to simulate missing traits");
        throw new RuntimeException("Not yet implemented");
    }

    private void computeOuterProducts() {

        // Make sure partials on tree are ready
        dataLikelihood.getLogLikelihood();

        if (likelihoodDelegate != outerProductDelegate) {
            simulateMissingTraits();
        }

        treeTraversalDelegate.updateAllNodes();
        treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

        final List<DataLikelihoodDelegate.BranchOperation> branchOperations = treeTraversalDelegate.getBranchOperations();
        final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = treeTraversalDelegate.getNodeOperations();

        final NodeRef root = dataLikelihood.getTree().getRoot();

        try {
            outerProductDelegate.setComputeWishartStatistics(true);
            outerProductDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());
            outerProductDelegate.setComputeWishartStatistics(false);

        } catch (DataLikelihoodDelegate.LikelihoodUnderflowException e) {
            throw new RuntimeException("Unhandled exception");
        }

        wishartStatistics = outerProductDelegate.getWishartStatistics();
    }

    @Override
    public MatrixParameterInterface getPrecisionParamter() {
        return likelihoodDelegate.getDiffusionModel().getPrecisionParameter();
    }

    @Override
    protected void storeState() {
        savedTraitDataKnown = traitDataKnown;
        savedOuterProductsKnown = outerProductsKnown;

        if (outerProductsKnown) {
            if (savedWishartStatistics == null) {
                savedWishartStatistics = wishartStatistics.clone();
            } else {
                wishartStatistics.copyTo(savedWishartStatistics);
            }
        }
    }

    @Override
    protected void restoreState() {
        traitDataKnown = savedTraitDataKnown;
        outerProductsKnown = savedOuterProductsKnown;

        if (outerProductsKnown) {
            WishartSufficientStatistics tmp = wishartStatistics;
            wishartStatistics = savedWishartStatistics;
            savedWishartStatistics = tmp;
        }
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        outerProductsKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        outerProductsKnown = false;
        // TODO If no partially missing traits and diffusion model hit, then no update necessary
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.hasId() ? xo.getId() : PARSER_NAME;
            String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

            final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

            if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("May not provide a sequence data likelihood in the precision Gibbs sampler");

            }

            final ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;

            return new WishartStatisticsWrapper(name,
                    traitName,
                    treeDataLikelihood, continuousData);
        }

        /**
         * @return an array of syntax rules required by this element.
         * Order is not important.
         */
        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return syntax;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return ConjugateWishartStatisticsProvider.class;
        }

        /**
         * @return Parser name, which is identical to name of xml element parsed by it.
         */
        @Override
        public String getParserName() {
            return PARSER_NAME;
        }

        private final XMLSyntaxRule[] syntax = new XMLSyntaxRule[] {
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TRAIT_NAME, true),
        };
    };


    private final ContinuousTraitDataModel continuousTraitDataModel = null;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final ContinuousDataLikelihoodDelegate outerProductDelegate;
    private final TreeDataLikelihood dataLikelihood;

    private final boolean hasPartiallyMissingTraits;

    private boolean traitDataKnown;
    private boolean outerProductsKnown;

    private boolean savedTraitDataKnown;
    private boolean savedOuterProductsKnown;

    private WishartSufficientStatistics wishartStatistics;
    private WishartSufficientStatistics savedWishartStatistics;

}
