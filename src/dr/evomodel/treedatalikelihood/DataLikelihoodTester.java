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

package dr.evomodel.treedatalikelihood;


import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodelxml.siteratemodel.GammaSiteModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
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
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DataLikelihoodTester {

    public static void main(String[] args) {

        SimpleAlignment alignment = createAlignment(sequences, Nucleotides.INSTANCE);

        TreeModel treeModel;

        try {
            treeModel = createSpecifiedTree("((human:0.1,chimp:0.1):0.1,gorilla:0.2)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        System.out.print("\nTest BeagleTreeLikelihood: ");

        //substitutionModel
        Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 1.0, 0, 100);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
//        double alpha = 0.5;
//        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("gammaModel", alpha, 4);
        GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteRateModel");
        siteRateModel.setSubstitutionModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.SUBSTITUTION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteRateModel.setMutationRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        BranchModel branchModel = new HomogeneousBranchModel(
                siteRateModel.getSubstitutionModel(),
                siteRateModel.getSubstitutionModel().getFrequencyModel());

        BranchRateModel branchRateModel = null;

        BeagleTreeLikelihood treeLikelihood = new BeagleTreeLikelihood(
                patterns,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                null,
                false, PartialsRescalingScheme.AUTO, true);

        double logLikelihood = treeLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("\nTest BeagleDataLikelihoodDelegate: ");

        BeagleDataLikelihoodDelegate dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                branchModel,
                patterns,
                siteRateModel, false);

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        MultiPartitionDataLikelihoodDelegate multiPartitionDataLikelihoodDelegate;

//        System.out.println("\nTest MultiPartitionDataLikelihoodDelegate 1 partition:");
//
//        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
//                treeModel,
//                Collections.singletonList((PatternList)patterns),
//                Collections.singletonList((BranchModel)branchModel),
//                Collections.singletonList((SiteRateModel)siteRateModel),
//                true);
//
//        treeDataLikelihood = new TreeDataLikelihood(
//                multiPartitionDataLikelihoodDelegate,
//                treeModel,
//                branchRateModel);
//
//        logLikelihood = treeDataLikelihood.getLogLikelihood();
//
//        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("\nTest MultiPartitionDataLikelihoodDelegate 2 partition: ");

        List<PatternList> patternLists = new ArrayList<>();
        patternLists.add(patterns);
        patternLists.add(patterns);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                Collections.singletonList(branchModel),
                Collections.singletonList((SiteRateModel)siteRateModel),
                true);

//        FrequencyModel f2 = new FrequencyModel(Nucleotides.INSTANCE, freqs);
//        HKY hky2 = new HKY(kappa, f2);
//
//        GammaSiteRateModel siteRateModel2 = new GammaSiteRateModel("siteRateModel");
//        siteRateModel2.setSubstitutionModel(hky2);
//        siteRateModel2.setMutationRateParameter(mu);
//
//        List<SiteRateModel> siteRateModels = new ArrayList<>();
//        siteRateModels.add(siteRateModel);
//        siteRateModels.add(siteRateModel2);
//
//        BranchModel branchModel2 = new HomogeneousBranchModel(
//                siteRateModel2.getSubstitutionModel(),
//                siteRateModel2.getSubstitutionModel().getFrequencyModel());
//
//        List<BranchModel> branchModels = new ArrayList<>();
//        branchModels.add(branchModel);
//        branchModels.add(branchModel2);
//
//        MultiPartitionDataLikelihoodDelegate multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
//                treeModel,
//                patternLists,
//                branchModels,
//                siteRateModels,
//                true);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);
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

}
