/*
 * EpochGradientRequiresAugmentedNodesTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.SubstitutionModelDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightGradientForDiscreteTrait;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;

/**
 * Branch-rate and node-height gradients on a convolving (epoch) branch model are only exact with augmented epoch
 * nodes. Without them the delegate must refuse to build instead of returning the gradient of the last epoch's Q.
 */
public class EpochGradientRequiresAugmentedNodesTest extends TraceCorrelationAssert {

    public EpochGradientRequiresAugmentedNodesTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testBranchRateGradientRefusedWithoutAugmentedNodes() {
        Fixture f = new Fixture(false);
        try {
            new BranchRateGradientForDiscreteTrait("test", f.likelihood, f.delegate, f.rates, false);
            fail("branch-rate gradient on an epoch model without augmented epoch nodes must be refused");
        } catch (IllegalArgumentException e) {
            assertEquals(SubstitutionModelDelegate.EPOCH_GRADIENT_MESSAGE, e.getMessage());
        }
    }

    public void testNodeHeightGradientRefusedWithoutAugmentedNodes() {
        Fixture f = new Fixture(false);
        try {
            new NodeHeightGradientForDiscreteTrait("test", f.likelihood, f.delegate, f.rates);
            fail("node-height gradient on an epoch model without augmented epoch nodes must be refused");
        } catch (IllegalArgumentException e) {
            assertEquals(SubstitutionModelDelegate.EPOCH_GRADIENT_MESSAGE, e.getMessage());
        }
    }

    public void testBranchRateGradientAllowedWithAugmentedNodes() {
        Fixture f = new Fixture(true);
        double[] gradient = new BranchRateGradientForDiscreteTrait("test", f.likelihood, f.delegate, f.rates, false)
                .getGradientLogDensity();
        assertEquals(f.rates.getDimension(), gradient.length);
    }

    private class Fixture {

        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final Parameter rates;

        Fixture(boolean augmented) {
            DefaultTreeModel tree = new DefaultTreeModel(ladderTree());

            List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();
            models.add(new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.1, 0.2, 0.3, 0.4})));
            models.add(new HKY(5.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.4, 0.3, 0.2, 0.1})));
            EpochBranchModel branchModel = new EpochBranchModel(tree, models, new Parameter.Default(new double[]{0.15}));

            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, 4);
            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

            rates = new Parameter.Default(tree.getNodeCount() - 1, 1.0);
            ArbitraryBranchRates branchRates = new ArbitraryBranchRates(tree, rates,
                    ArbitraryBranchRates.make(false, false, false), false);

            delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel, siteRateModel,
                    false, false, PartialsRescalingScheme.NONE, false,
                    new PreOrderSettings(true, true, false, false, false, false, augmented));
            likelihood = new TreeDataLikelihood(delegate, tree, branchRates);
        }

        private SimpleTree ladderTree() {
            SimpleNode current = new SimpleNode();
            current.setTaxon(taxa[0]);
            current.setHeight(0.0);
            for (int i = 1; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                SimpleNode parent = new SimpleNode();
                parent.setHeight(0.05 * i);
                parent.addChild(current);
                parent.addChild(tip);
                current = parent;
            }
            return new SimpleTree(current);
        }
    }
}
