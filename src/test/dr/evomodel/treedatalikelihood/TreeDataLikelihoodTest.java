/*
 * LikelihoodTest.java
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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.GTRParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;

import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;


/**
 * @author Andrew Rambuat
 * convert testLikelihood.xml in the folder /example
 */

public class TreeDataLikelihoodTest extends TraceCorrelationAssert {

    private TreeModel treeModel;
    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public TreeDataLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        treeModel = createPrimateTreeModel ();
    }

    public void testNewickTree() {
        System.out.println("\nTest Simple Node to convert Newick Tree:");
        String expectedNewickTree = "((((human:0.024003,(chimp:0.010772,bonobo:0.010772):0.013231):0.012035," +
                "gorilla:0.036038):0.033087,orangutan:0.069125):0.030457,siamang:0.099582);";

        assertEquals("Fail to covert the correct tree !!!", expectedNewickTree, TreeUtils.newick(treeModel, 6));
    }

    private TreeDataLikelihood getTreeDataLikelihood(SubstitutionModel substitutionModel, GammaSiteRateModel siteRateModel) {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        DataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                new HomogeneousBranchModel(substitutionModel),
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.DEFAULT,
                false,
                new PreOrderSettings(false, false, false, true)
        );

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                new DefaultBranchRateModel());
        return treeDataLikelihood;
    }


    public void testLikelihoodJC69() {
        System.out.println("\nTest Likelihood using JC69:");
        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0);

        FrequencyModel frequencyModel = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, frequencyModel);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma");

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodJC69", format.format(-1992.20564), format.format(treeDataLikelihood.getLogLikelihood()));
    }


    public void testLikelihoodK80() {
        System.out.println("\nTest Likelihood using K80:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25}));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma");

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodK80", format.format(-1868.89782), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85() {
        System.out.println("\nTest Likelihood using HKY85:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma");

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodHKY85", format.format(-1839.84514), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85G() {
        System.out.println("\nTest Likelihood using HKY85G:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                new Parameter.Default(0.5),
                4, null);

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodHKY85G", format.format(-1816.82611), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85I() {
        System.out.println("\nTest Likelihood using HKY85I:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                null,
                1,
                new Parameter.Default(0.75));

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodHKY85I", format.format(-1822.37478), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85GI() {
        System.out.println("\nTest Likelihood using HKY85GI:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 8.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                new Parameter.Default(0.5),
                4,
                new Parameter.Default(0.75));

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(hky, siteRateModel);

        assertEquals("treeLikelihoodHKY85GI", format.format(-1815.02176), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTR() {
        System.out.println("\nTest Likelihood using GTR:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        GTR gtr = new GTR(
                new Parameter.Default(GTRParser.A_TO_C, 1.0),
                new Parameter.Default(GTRParser.A_TO_G, 1.0),
                new Parameter.Default(GTRParser.A_TO_T, 1.0),
                new Parameter.Default(GTRParser.C_TO_G, 1.0),
                new Parameter.Default(GTRParser.C_TO_T, 1.0),
                new Parameter.Default(GTRParser.G_TO_T, 1.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                null, 1,
                null);

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(gtr, siteRateModel);

        assertEquals("treeLikelihoodGTR", format.format(-1969.14584), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRI() {
        System.out.println("\nTest Likelihood using GTRI:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        GTR gtr = new GTR(
                new Parameter.Default(GTRParser.A_TO_C, 1.0),
                new Parameter.Default(GTRParser.A_TO_G, 1.0),
                new Parameter.Default(GTRParser.A_TO_T, 1.0),
                new Parameter.Default(GTRParser.C_TO_G, 1.0),
                new Parameter.Default(GTRParser.C_TO_T, 1.0),
                new Parameter.Default(GTRParser.G_TO_T, 1.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                null, 1,
                new Parameter.Default(0.5));

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(gtr, siteRateModel);

        assertEquals("treeLikelihoodGTRI", format.format(-1948.84175), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRG() {
        System.out.println("\nTest Likelihood using GTRG:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        GTR gtr = new GTR(
                new Parameter.Default(GTRParser.A_TO_C, 1.0),
                new Parameter.Default(GTRParser.A_TO_G, 1.0),
                new Parameter.Default(GTRParser.A_TO_T, 1.0),
                new Parameter.Default(GTRParser.C_TO_G, 1.0),
                new Parameter.Default(GTRParser.C_TO_T, 1.0),
                new Parameter.Default(GTRParser.G_TO_T, 1.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                new Parameter.Default(0.5), 4,
                null);

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(gtr, siteRateModel);

        assertEquals("treeLikelihoodGTRG", format.format(-1949.03601), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRGI() {
        System.out.println("\nTest Likelihood using GTRGI:");

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(alignment.getStateFrequencies()));
        GTR gtr = new GTR(
                new Parameter.Default(GTRParser.A_TO_C, 1.0),
                new Parameter.Default(GTRParser.A_TO_G, 1.0),
                new Parameter.Default(GTRParser.A_TO_T, 1.0),
                new Parameter.Default(GTRParser.C_TO_G, 1.0),
                new Parameter.Default(GTRParser.C_TO_T, 1.0),
                new Parameter.Default(GTRParser.G_TO_T, 1.0), f);

        //siteModel
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gamma",
                null, 1.0,
                new Parameter.Default(0.5), 4,
                new Parameter.Default(0.5));

        //treeLikelihood
        TreeDataLikelihood treeDataLikelihood = getTreeDataLikelihood(gtr, siteRateModel);

        assertEquals("treeLikelihoodGTRGI", format.format(-1951.62188), format.format(treeDataLikelihood.getLogLikelihood()));
    }

    public static Test suite() {
        return new TestSuite(TreeDataLikelihoodTest.class);
    }
}
