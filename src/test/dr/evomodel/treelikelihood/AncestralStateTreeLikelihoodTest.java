/*
 * AncestralStateTreeLikelihoodTest.java
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

package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.HomogeneousRateDelegate;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.Map;
import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class AncestralStateTreeLikelihoodTest extends TestCase {

    private FlexibleTree tree;

    public void setUp() throws Exception {
        super.setUp();

        MathUtils.setSeed(666);

        NewickImporter importer = new NewickImporter("(0:2.0,(1:1.0,2:1.0):1.0);");
        tree = (FlexibleTree) importer.importTree(null);
    }

    // transition prob of JC69
    private double t(boolean same, double time) {
        if (same) {
            return 0.25 + 0.75 * Math.exp(-4.0 / 3.0 * time);
        } else {
            return 0.25 - 0.25 * Math.exp(-4.0 / 3.0 * time);
        }
    }

    public void testJointLikelihood() {

        TreeModel treeModel = new DefaultTreeModel("treeModel", tree);

        Sequence[] sequence = new Sequence[3];

        sequence[0] = new Sequence(new Taxon("0"), "A");
        sequence[1] = new Sequence(new Taxon("1"), "C");
        sequence[2] = new Sequence(new Taxon("2"), "C");

        Taxa taxa = new Taxa();
        for (Sequence s : sequence) {
            taxa.addTaxon(s.getTaxon());
        }

        SimpleAlignment alignment = new SimpleAlignment();
        for (Sequence s : sequence) {
            alignment.addSequence(s);
        }

        Parameter mu = new Parameter.Default(1, 1.0);

        Parameter kappa = new Parameter.Default(1, 1.0);
        double[] pi = {0.25, 0.25, 0.25, 0.25};

        Parameter freqs = new Parameter.Default(pi);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        AncestralStateBeagleTreeLikelihood treeLikelihood = new AncestralStateBeagleTreeLikelihood(
                alignment,
                treeModel,
                new HomogeneousBranchModel(hky),
                new DiscretizedSiteRateModel(""),
                new StrictClockBranchRates(mu),
                null,
                false, PartialsRescalingScheme.DEFAULT,
                false,
                null,
                Nucleotides.INSTANCE,
                "state",
                false,
                true, // useMap = true
                false);

        double logLike = treeLikelihood.getLogLikelihood();

        StringBuffer buffer = new StringBuffer();

        TreeUtils.newick(treeModel, treeModel.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                null, null, new TreeTraitProvider[] { treeLikelihood }, null, buffer);


        System.out.println(buffer);

        System.out.println("t_CA(2) = " + t(false, 2.0));
        System.out.println("t_CC(1) = " + t(true, 1.0));

        double trueValue = 0.25 * t(false, 2.0) * Math.pow(t(true, 1.0), 3.0);

        assertEquals(logLike, Math.log(trueValue), 1e-6);
    }
}
