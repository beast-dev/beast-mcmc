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

package test.dr.evolution;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
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

public class SitePatternsTest extends TraceCorrelationAssert {

    public SitePatternsTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testAllSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.UNCOMPRESSED);
        System.out.println("Uncompressed site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 768);
    }

    public void testUniqueSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment);
        System.out.println("Unique site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 69);
    }

    public void testAmbConstSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.AMBIGUOUS_CONSTANT);
        System.out.println("Unique + ambiguous constant patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 68);
    }

    public void testConstantOnlySitePatterns() {
        SitePatterns uncompressed = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.UNCOMPRESSED);
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.CONSTANT_ONLY);
        int stateCount = patterns.getStateCount();

        // the first stateCount patterns are the constant ones, weighted by their number of sites
        double[] constantWeights = new double[stateCount];
        int variableSites = 0;
        for (int site = 0; site < uncompressed.getPatternCount(); site++) {
            int[] pattern = uncompressed.getPattern(site);
            int constantState = constantState(pattern, stateCount);
            if (constantState >= 0) {
                constantWeights[constantState] += 1;
                assertEquals(constantState, patterns.getPatternIndex(site));
            } else {
                int index = patterns.getPatternIndex(site);
                assertTrue(index >= stateCount);
                assertTrue(java.util.Arrays.equals(pattern, patterns.getPattern(index)));
                variableSites++;
            }
        }
        assertEquals(stateCount + variableSites, patterns.getPatternCount());
        for (int i = 0; i < stateCount; i++) {
            assertEquals(i, constantState(patterns.getPattern(i), stateCount));
            assertEquals(constantWeights[i], patterns.getPatternWeight(i), 0.0);
        }
        for (int i = stateCount; i < patterns.getPatternCount(); i++) {
            assertEquals(1.0, patterns.getPatternWeight(i), 0.0);
        }
    }

    public void testAmbiguousConstantOnlySitePatterns() {
        // columns: constant A; A with one N; A with one R; A with one Y; A with three Ns; constant C; variable
        String[] columns = {"AAAAAA", "AANAAA", "AARAAA", "AAYAAA", "ANNNAA", "CCCCCC", "ACAAAA"};
        String[] taxa = {"t1", "t2", "t3", "t4", "t5", "t6"};
        String[] sequences = new String[taxa.length];
        for (int t = 0; t < taxa.length; t++) {
            StringBuilder sb = new StringBuilder();
            for (String column : columns) {
                sb.append(column.charAt(t));
            }
            sequences[t] = sb.toString();
        }
        createAlignment(new Object[][]{taxa, sequences}, Nucleotides.INSTANCE);

        SitePatterns constantOnly = new SitePatterns(alignment, null, 0, -1, 1, false, SitePatterns.CompressionType.CONSTANT_ONLY);
        assertEquals(4 + 5, constantOnly.getPatternCount());
        assertEquals(1.0, constantOnly.getPatternWeight(Nucleotides.A_STATE), 0.0);
        assertEquals(1.0, constantOnly.getPatternWeight(Nucleotides.C_STATE), 0.0);

        // with the default threshold (25%), one N or a compatible R still counts as constant, but not a Y or three Ns
        SitePatterns ambiguousConstantOnly = new SitePatterns(alignment, null, 0, -1, 1, false, SitePatterns.CompressionType.AMBIGUOUS_CONSTANT_ONLY);
        assertEquals(4 + 3, ambiguousConstantOnly.getPatternCount());
        assertEquals(3.0, ambiguousConstantOnly.getPatternWeight(Nucleotides.A_STATE), 0.0);
        assertEquals(1.0, ambiguousConstantOnly.getPatternWeight(Nucleotides.C_STATE), 0.0);
        assertEquals(Nucleotides.A_STATE, ambiguousConstantOnly.getPatternIndex(1));
        assertEquals(Nucleotides.A_STATE, ambiguousConstantOnly.getPatternIndex(2));
        assertTrue(ambiguousConstantOnly.getPatternIndex(3) >= 4);
        assertTrue(ambiguousConstantOnly.getPatternIndex(4) >= 4);

        // the existing AMBIGUOUS_CONSTANT compression must apply the same threshold: the site with three Ns
        // (50% ambiguous) is not merged into the constant A pattern
        SitePatterns ambiguousConstant = new SitePatterns(alignment, null, 0, -1, 1, false, SitePatterns.CompressionType.AMBIGUOUS_CONSTANT);
        assertEquals(3.0, ambiguousConstant.getPatternWeight(Nucleotides.A_STATE), 0.0);
    }

    private static int constantState(int[] pattern, int stateCount) {
        for (int state : pattern) {
            if (state != pattern[0] || state >= stateCount) {
                return -1;
            }
        }
        return pattern[0];
    }

    public void testCodonSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 2, -1, 3);
        System.out.println("2nd position unique site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 37);
    }


    public static Test suite() {
        return new TestSuite(SitePatternsTest.class);
    }
}
