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

package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.TreeUtils;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.GTR;
import dr.oldevomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.TreeLikelihood;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.oldevomodelxml.substmodel.GTRParser;
import dr.oldevomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;


/**
 * @author Walter Xie
 * convert testLikelihood.xml in the folder /example
 */

public class LikelihoodTest extends TraceCorrelationAssert {

    private TreeModel treeModel;
    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public LikelihoodTest(String name) {
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



    public void testLikelihoodJC69() {
        System.out.println("\nTest Likelihood using JC69:");
        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodJC69", format.format(-1992.20564), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodK80() {
        System.out.println("\nTest Likelihood using K80:");
        // Sub model
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 27.402591, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodK80", format.format(-1856.30305), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85() {
        System.out.println("\nTest Likelihood using HKY85:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 29.739445, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);          

        assertEquals("treeLikelihoodHKY85", format.format(-1825.21317), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85G() {
        System.out.println("\nTest Likelihood using HKY85G:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 38.829740, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModelParser.GAMMA_SHAPE, 0.137064, 0, 1000.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, shape, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85G", format.format(-1789.75936), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85I() {
        System.out.println("\nTest Likelihood using HKY85I:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 38.564672, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter invar = new Parameter.Default(GammaSiteModelParser.PROPORTION_INVARIANT, 0.701211, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, null, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85I", format.format(-1789.91240), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodHKY85GI() {
        System.out.println("\nTest Likelihood using HKY85GI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 39.464538, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModelParser.GAMMA_SHAPE, 0.587649, 0, 1000.0);
        Parameter invar = new Parameter.Default(GammaSiteModelParser.PROPORTION_INVARIANT, 0.486548, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(hky, mu, shape, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodHKY85GI", format.format(-1789.63923), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTR() {
        System.out.println("\nTest Likelihood using GTR:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTRParser.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTRParser.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTRParser.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTRParser.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTRParser.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTRParser.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, null, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTR", format.format(-1969.14584), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRI() {
        System.out.println("\nTest Likelihood using GTRI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTRParser.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTRParser.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTRParser.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTRParser.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTRParser.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTRParser.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter invar = new Parameter.Default(GammaSiteModelParser.PROPORTION_INVARIANT, 0.5, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, null, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRI", format.format(-1948.84175), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRG() {
        System.out.println("\nTest Likelihood using GTRG:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTRParser.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTRParser.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTRParser.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTRParser.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTRParser.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTRParser.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModelParser.GAMMA_SHAPE, 0.5, 0, 100.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, shape, 4, null);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRG", format.format(-1949.03601), format.format(treeLikelihood.getLogLikelihood()));
    }

    public void testLikelihoodGTRGI() {
        System.out.println("\nTest Likelihood using GTRGI:");
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);

        Variable<Double> rateACValue = new Parameter.Default(GTRParser.A_TO_C, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateAGValue = new Parameter.Default(GTRParser.A_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateATValue = new Parameter.Default(GTRParser.A_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCGValue = new Parameter.Default(GTRParser.C_TO_G, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateCTValue = new Parameter.Default(GTRParser.C_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        Variable<Double> rateGTValue = new Parameter.Default(GTRParser.G_TO_T, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        GTR gtr = new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, f);

        //siteModel
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        Parameter shape = new Parameter.Default(GammaSiteModelParser.GAMMA_SHAPE, 0.5, 0, 100.0);
        Parameter invar = new Parameter.Default(GammaSiteModelParser.PROPORTION_INVARIANT, 0.5, 0, 1.0);

        GammaSiteModel siteModel = new GammaSiteModel(gtr, mu, shape, 4, invar);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        assertEquals("treeLikelihoodGTRGI", format.format(-1947.58294), format.format(treeLikelihood.getLogLikelihood()));
    }


    public static Test suite() {
        return new TestSuite(LikelihoodTest.class);
    }
}
