/*
 * PredictiveBeagleSequenceSimulatorGeneratorTest.java
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

package test.dr.evomodel.treedatalikelihood.pps;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.pps.PredictiveBeagleSequenceSimulatorGenerator;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class PredictiveBeagleSequenceSimulatorGeneratorTest extends TestCase {

    private SimpleAlignment alignment;
    private TreeModel treeModel;
    private SitePatterns patterns;
    private HomogeneousBranchModel branchModel;
    private FrequencyModel frequencyModel;
    private GammaSiteRateModel singleCategorySiteRateModel;

    public PredictiveBeagleSequenceSimulatorGeneratorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        alignment = createAlignment();

        NewickImporter importer = new NewickImporter("(A:0.1,(B:0.2,C:0.15):0.1);");
        treeModel = new DefaultTreeModel((FlexibleTree) importer.importTree(null));

        patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        frequencyModel = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(alignment.getStateFrequencies()));
        HKY hky = new HKY(new Parameter.Default(HKYParser.KAPPA, 2.0), frequencyModel);
        branchModel = new HomogeneousBranchModel(hky);

        singleCategorySiteRateModel = new GammaSiteRateModel("gamma");
    }

    private static SimpleAlignment createAlignment() {
        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);

        addSequence(alignment, "A", "ACGT");
        addSequence(alignment, "B", "ACGA");
        addSequence(alignment, "C", "ACGC");

        return alignment;
    }

    private static void addSequence(SimpleAlignment alignment, String taxonId, String sequenceString) {
        Sequence sequence = new Sequence(sequenceString);
        sequence.setTaxon(new Taxon(taxonId));
        sequence.setDataType(Nucleotides.INSTANCE);
        alignment.addSequence(sequence);
    }

    private TreeDataLikelihood createTreeDataLikelihood(PatternList patternList, GammaSiteRateModel siteRateModel) {
        DataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patternList,
                branchModel,
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.DEFAULT,
                false,
                new PreOrderSettings(false, false, false, true));

        return new TreeDataLikelihood(delegate, treeModel, new DefaultBranchRateModel());
    }

    public void testSimulateProducesExpectedTaxaAndSiteCount() {
        TreeDataLikelihood treeDataLikelihood = createTreeDataLikelihood(patterns, singleCategorySiteRateModel);

        PredictiveBeagleSequenceSimulatorGenerator generator =
                new PredictiveBeagleSequenceSimulatorGenerator(patterns, treeDataLikelihood, false);

        SimpleAlignment simulated = generator.simulate();

        assertEquals(alignment.getTaxonCount(), simulated.getTaxonCount());
        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            assertTrue("simulated alignment is missing taxon " + alignment.getTaxonId(i),
                    simulated.getTaxonIndex(alignment.getTaxonId(i)) != -1);
        }
        for (int i = 0; i < simulated.getSequenceCount(); i++) {
            assertEquals(alignment.getSiteCount(), simulated.getSequence(i).getLength());
        }

        assertEquals(Nucleotides.INSTANCE, generator.getDataType());
        assertEquals(patterns.getId(), generator.getParentId());
    }

    public void testSimulateThrowsWhenPatternListHasTaxaNotInTheTree() {
        SimpleAlignment biggerAlignment = new SimpleAlignment();
        biggerAlignment.setDataType(Nucleotides.INSTANCE);
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            biggerAlignment.addSequence(alignment.getSequence(i));
        }
        addSequence(biggerAlignment, "extraTaxonNotInTree", alignment.getSequence(0).getSequenceString());

        SitePatterns biggerPatterns = new SitePatterns(biggerAlignment, null, 0, -1, 1, true);

        TreeDataLikelihood treeDataLikelihood = createTreeDataLikelihood(biggerPatterns, singleCategorySiteRateModel);

        PredictiveBeagleSequenceSimulatorGenerator generator =
                new PredictiveBeagleSequenceSimulatorGenerator(biggerPatterns, treeDataLikelihood, false);

        try {
            generator.simulate();
            fail("Expected a RuntimeException because the tree has fewer taxa than the expected pattern list");
        } catch (RuntimeException e) {
            // expected
        }
    }

    public static Test suite() {
        return new TestSuite(PredictiveBeagleSequenceSimulatorGeneratorTest.class);
    }
}
