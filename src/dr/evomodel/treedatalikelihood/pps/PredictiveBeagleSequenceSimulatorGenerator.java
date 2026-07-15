/*
 * PredictiveBeagleSequenceSimulatorGenerator.java
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

package dr.evomodel.treedatalikelihood.pps;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

import java.util.ArrayList;

/**
 * Hijack piBUSS machinery for PPS
 */
public class PredictiveBeagleSequenceSimulatorGenerator implements PredictiveDataGenerator {

    private final TreeModel treeModel;
    private final BranchModel branchModel;
    private final GammaSiteRateModel siteRateModel;
    private final BranchRateModel branchRateModel;
    private final FrequencyModel freqModel;
    private final PatternList patternList;
    private final DataType dataType;
    private final boolean outputAncestralSequences;

    public PredictiveBeagleSequenceSimulatorGenerator(PatternList realDataPatterns, TreeDataLikelihood treeDataLikelihood,
                                                       boolean outputAncestralSequences) {

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new RuntimeException("predictiveBeagleSequenceSimulator requires a BeagleDataLikelihoodDelegate");
        }
        BeagleDataLikelihoodDelegate beagleDelegate = (BeagleDataLikelihoodDelegate) delegate;

        PatternList delegatePatterns = beagleDelegate.getPatternList();
        if (realDataPatterns.getPatternCount() != delegatePatterns.getPatternCount() ||
                !realDataPatterns.asList().containsAll(delegatePatterns.asList()) ||
                !delegatePatterns.asList().containsAll(realDataPatterns.asList())) {
            throw new RuntimeException("predictiveBeagleSequenceSimulator's supplied patterns '" +
                    realDataPatterns.getId() + "' is not consistent with the patterns used by the " +
                    "treeDataLikelihood's BeagleDataLikelihoodDelegate (different taxa or pattern count)");
        }

        SiteRateModel resolvedSiteRateModel = beagleDelegate.getSiteRateModel();
        if (!(resolvedSiteRateModel instanceof GammaSiteRateModel) || resolvedSiteRateModel.getCategoryCount() > 1) {
            throw new RuntimeException("predictiveBeagleSequenceSimulator requires a single-category " +
                    "GammaSiteRateModel (among-site rate variation is not supported), but found " +
                    resolvedSiteRateModel.getClass().getName() + " with " +
                    resolvedSiteRateModel.getCategoryCount() + " categories");
        }
        this.siteRateModel = (GammaSiteRateModel) resolvedSiteRateModel;

        this.patternList = realDataPatterns;
        this.dataType = patternList.getDataType();

        this.branchModel = beagleDelegate.getBranchModel();
        this.freqModel = branchModel.getRootFrequencyModel();
        this.branchRateModel = treeDataLikelihood.getBranchRateModel();

        Tree tree = treeDataLikelihood.getTree();
        if (!(tree instanceof TreeModel)) {
            throw new RuntimeException("predictiveBeagleSequenceSimulator requires the treeDataLikelihood's " +
                    "tree to be a TreeModel, but found " + tree.getClass().getName());
        }
        this.treeModel = (TreeModel) tree;
        this.outputAncestralSequences = outputAncestralSequences;
    }

    @Override
    public SimpleAlignment simulate() {

        int siteCount = getExpectedSiteCount();

        Partition partition = new Partition(treeModel, branchModel, siteRateModel, branchRateModel, freqModel,
                0, siteCount - 1, 1, dataType);

        ArrayList<Partition> partitions = new ArrayList<Partition>();
        partitions.add(partition);
        BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(partitions);
        SimpleAlignment alignment = simulator.simulate(false, outputAncestralSequences);

        validate(alignment, siteCount);

        return alignment;
    }

    private int getExpectedSiteCount() {
        double weightSum = 0.0;
        for (double weight : patternList.getPatternWeights()) {
            weightSum += weight;
        }
        return (int) Math.round(weightSum);
    }

    private void validate(SimpleAlignment alignment, int expectedSiteCount) {

        TaxonList expectedTaxa = getTaxa();

        if (alignment.getTaxonCount() != expectedTaxa.getTaxonCount()) {
            throw new RuntimeException("Simulated alignment has " + alignment.getTaxonCount() +
                    " taxa but expected " + expectedTaxa.getTaxonCount());
        }

        for (int i = 0; i < expectedTaxa.getTaxonCount(); i++) {
            Taxon taxon = expectedTaxa.getTaxon(i);
            if (alignment.getTaxonIndex(taxon.getId()) == -1) {
                throw new RuntimeException("Simulated alignment is missing expected taxon " + taxon.getId());
            }
        }

        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            int length = alignment.getSequence(i).getLength();
            if (length != expectedSiteCount) {
                throw new RuntimeException("Simulated sequence for " + alignment.getTaxon(i).getId() +
                        " has length " + length + " but expected " + expectedSiteCount);
            }
        }
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    @Override
    public String getParentId() {
        return patternList.getId();
    }

    @Override
    public TaxonList getTaxa() {
        return patternList;
    }

    public FrequencyModel getFrequencyModel() {
        return freqModel;
    }
}
