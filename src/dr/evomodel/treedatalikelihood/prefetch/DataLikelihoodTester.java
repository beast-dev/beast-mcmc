/*
 * MultiPartitionDataLikelihoodDelegate.java
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

package dr.evomodel.treedatalikelihood.prefetch;


import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.operators.PrefetchSubtreeLeapOperator;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.siteratemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.math.MathUtils;

import java.util.Collections;
import java.util.logging.Logger;

public class DataLikelihoodTester {

    public static void main(String[] args) {

        MathUtils.setSeed(666);
        
        // turn off logging to avoid screen noise...
        Logger logger = Logger.getLogger("dr");
        logger.setUseParentHandlers(false);

        SimpleAlignment alignment = createAlignment(sequences, Nucleotides.INSTANCE);

        TreeModel treeModel;

        try {
            treeModel = createSpecifiedTree("((human:0.1,chimp:0.1):0.1,gorilla:0.2)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        //substitutionModel
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        double alpha = 0.5;
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gammaModel", alpha, 4);
//        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteRateModel");
        siteRateModel.setSubstitutionModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.SUBSTITUTION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteRateModel.setRelativeRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        BranchModel branchModel = new HomogeneousBranchModel(
                siteRateModel.getSubstitutionModel(),
                siteRateModel.getSubstitutionModel().getFrequencyModel());

        BranchRateModel branchRateModel = new DefaultBranchRateModel();

        BeagleDataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                branchModel,
                siteRateModel, false,
                PartialsRescalingScheme.NONE,
                false);

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        double logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println();
        System.out.println("BeagleDataLikelihoodDelegate:         logLikelihood = " + logLikelihood);

        MultiPartitionDataLikelihoodDelegate multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)patterns),
                Collections.singletonList(branchModel),
                Collections.singletonList((SiteRateModel)siteRateModel),
                false,
                PartialsRescalingScheme.NONE,
                false);

        TreeDataLikelihood treeDataLikelihood2 = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood2.getLogLikelihood();

        System.out.println("MultipartitionDataLikelihoodDelegate: logLikelihood = " + logLikelihood);

        PrefetchDataLikelihoodDelegate prefetchDataLikelihoodDelegate = new PrefetchDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)patterns),
                Collections.singletonList(branchModel),
                Collections.singletonList((SiteRateModel)siteRateModel),
                false,
                PartialsRescalingScheme.NONE,
                false,
                2);

        PrefetchTreeDataLikelihood prefetchTreeDataLikelihood = new PrefetchTreeDataLikelihood(
                prefetchDataLikelihoodDelegate,
                treeModel,
                branchRateModel);
        
        logLikelihood = prefetchTreeDataLikelihood.getLogLikelihood();

        System.out.println("PrefetchDataLikelihoodDelegate:       logLikelihood = " + logLikelihood);
        System.out.println();

//        treeDataLikelihood2.storeModelState();
        System.out.println("Store state");
        prefetchTreeDataLikelihood.storeModelState();

        PrefetchSubtreeLeapOperator op = new PrefetchSubtreeLeapOperator(prefetchTreeDataLikelihood, treeModel, 1.0, 1.0, 0.23, CoercionMode.COERCION_OFF);

        System.out.println("Operate - prefetch operations 1 & 2");
        op.operate();

//        System.out.println("Operation 1, MultiPartitionDataLikelihoodDelegate: lnL = " + treeDataLikelihood2.getLogLikelihood());
        System.out.println("Operation 1, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());

        System.out.println("Reject operation 1");
        op.reject();

//        treeDataLikelihood2.restoreModelState();
        System.out.println("Restore state (ignored)");
        prefetchTreeDataLikelihood.restoreModelState();

//        treeDataLikelihood2.storeModelState();
        System.out.println("Store state (ignored)");
        prefetchTreeDataLikelihood.storeModelState();

        System.out.println("Operate - switch to operation 2");
        op.operate();

//        System.out.println("Operation 2, MultiPartitionDataLikelihoodDelegate: lnL = " + treeDataLikelihood2.getLogLikelihood());
        System.out.println("Operation 2, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());

        System.out.println("Accept operation 2");
        op.accept(0.1);

        System.out.println("Accept state");
        prefetchTreeDataLikelihood.acceptModelState();

        System.out.println();

        //        treeDataLikelihood2.storeModelState();
        System.out.println("Store state");
        prefetchTreeDataLikelihood.storeModelState();

        System.out.println("Operate - prefetch operations 3 & 4");
        op.operate();

//        System.out.println("Operation 2, MultiPartitionDataLikelihoodDelegate: lnL = " + treeDataLikelihood2.getLogLikelihood());
        System.out.println("Operation 3, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());

        System.out.println("Accept operation 3");
        op.accept(0.1);

        System.out.println("Accept state");
        prefetchTreeDataLikelihood.acceptModelState();

        System.out.println();

        //        treeDataLikelihood2.storeModelState();
        System.out.println("Store state");
        prefetchTreeDataLikelihood.storeModelState();

        System.out.println("Operate - prefetch operations 5 & 6");
        op.operate();

//        System.out.println("Operation 2, MultiPartitionDataLikelihoodDelegate: lnL = " + treeDataLikelihood2.getLogLikelihood());
        System.out.println("Operation 5, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());

        System.out.println("Reject operation 5");
        op.reject();

//        treeDataLikelihood2.restoreModelState();
        System.out.println("Restore state (ignored)");
        prefetchTreeDataLikelihood.restoreModelState();

//        treeDataLikelihood2.storeModelState();
        System.out.println("Store state (ignored)");
        prefetchTreeDataLikelihood.storeModelState();

        System.out.println("Operate - switch to operation 6");
        op.operate();

//        System.out.println("Operation 2, MultiPartitionDataLikelihoodDelegate: lnL = " + treeDataLikelihood2.getLogLikelihood());
        System.out.println("Operation 6, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());

        System.out.println("Reject operation 6 (both prefetched 5 & 6 rejected)");
        op.reject();

        System.out.println();

//        treeDataLikelihood2.restoreModelState();
        System.out.println("Restore state");
        prefetchTreeDataLikelihood.restoreModelState();

        prefetchTreeDataLikelihood.suspendPrefetch();

        System.out.println("Restored lnL, PrefetchTreeDataLikelihood:           lnL = " + prefetchTreeDataLikelihood.getLogLikelihood());
        System.out.println();

    }

    private static SimpleAlignment createAlignment(Object[][] taxa_sequence, DataType dataType) {

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
//        alignment.setDataType(Nucleotides.INSTANCE);

        Taxon[] taxa = new Taxon[taxa_sequence[0].length]; // 6, 17
        System.out.println("Taxon len = " + taxa_sequence[0].length);
        System.out.println("Alignment len = " + taxa_sequence[1].length);
        if (taxa_sequence.length > 2) System.out.println("Date len = " + taxa_sequence[2].length);

        for (int i=0; i < taxa_sequence[0].length; i++) {
            taxa[i] = new Taxon(taxa_sequence[0][i].toString());

            if (taxa_sequence.length > 2) {
                Date date = new Date((Double) taxa_sequence[2][i], Units.Type.YEARS, (Boolean) taxa_sequence[3][0]);
                taxa[i].setDate(date);
            }

            //taxonList.addTaxon(taxon);
            Sequence sequence = new Sequence(taxa_sequence[1][i].toString());
            sequence.setTaxon(taxa[i]);
            sequence.setDataType(dataType);

            alignment.addSequence(sequence);
        }

        return alignment;

    }

    private static TreeModel createSpecifiedTree(String t) throws Exception {
        NewickImporter importer = new NewickImporter(t);
        Tree tree = importer.importTree(null);

        return new TreeModel(tree);//treeModel
    }

    static private String sequences[][] = {
            {"human", "chimp", "gorilla"},
            {
                    "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA",
                    "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA",
                    "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA"
            }
    };

    static private String moreSequences[][] = {
            {"human", "chimp", "gorilla"},
            {
                    "AGGGATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAAATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA",
                    "AGCGATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAAATTTCTGCTAGGTCTATACGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA",
                    "AGGTATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAAATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA"
            }
    };

}
