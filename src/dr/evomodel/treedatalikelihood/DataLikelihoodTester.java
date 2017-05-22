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


import dr.evolution.alignment.Patterns;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
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
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.util.MessageLogHandler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.*;

@Deprecated
public class DataLikelihoodTester {

    public static void main(String[] args) {

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

        System.out.print("\nTest BeagleTreeLikelihood (kappa = 1): ");

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

        FrequencyModel f2 = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        Parameter kappa2 = new Parameter.Default(HKYParser.KAPPA, 10.0, 0, 100);
        HKY hky2 = new HKY(kappa2, f2)
                ;

        GammaSiteRateModel siteRateModel2 = new GammaSiteRateModel("gammaModel", alpha, 4);
        siteRateModel2.setSubstitutionModel(hky2);
        siteRateModel2.setRelativeRateParameter(mu);

        //treeLikelihood
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

        BranchModel branchModel = new HomogeneousBranchModel(
                siteRateModel.getSubstitutionModel(),
                siteRateModel.getSubstitutionModel().getFrequencyModel());

        BranchModel branchModel2 = new HomogeneousBranchModel(
                siteRateModel2.getSubstitutionModel(),
                siteRateModel2.getSubstitutionModel().getFrequencyModel());

        BranchRateModel branchRateModel = new DefaultBranchRateModel();

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

        System.out.print("\nTest BeagleDataLikelihoodDelegate (kappa = 1): ");

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

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        hky.setKappa(5.0);
        System.out.print("\nTest BeagleDataLikelihoodDelegate (kappa = 5): ");
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("\nTest BeagleDataLikelihoodDelegate (kappa = 10): ");

        dataLikelihoodDelegate = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                branchModel2,
                siteRateModel2, false,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                dataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);
        hky2.setKappa(11.0);
        System.out.print("\nTest BeagleDataLikelihoodDelegate (kappa = 11): ");
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        hky.setKappa(1.0);
        hky2.setKappa(10.0);

        MultiPartitionDataLikelihoodDelegate multiPartitionDataLikelihoodDelegate;

        System.out.print("\nTest MultiPartitionDataLikelihoodDelegate 1 partition (kappa = 1):");

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)patterns),
                Collections.singletonList((BranchModel)branchModel),
                Collections.singletonList((SiteRateModel)siteRateModel),
                true,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        hky.setKappa(5.0);
        System.out.print("\nTest MultiPartitionDataLikelihoodDelegate 1 partition (kappa = 5):");
        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);
        hky.setKappa(1.0);

        System.out.print("\nTest MultiPartitionDataLikelihoodDelegate 1 partition (kappa = 10):");

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)patterns),
                Collections.singletonList((BranchModel)branchModel2),
                Collections.singletonList((SiteRateModel)siteRateModel2),
                true,
                PartialsRescalingScheme.NONE,
                false);


        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("\nTest MultiPartitionDataLikelihoodDelegate 2 partitions (kappa = 1, 10): ");

        List<PatternList> patternLists = new ArrayList<PatternList>();
        patternLists.add(patterns);
        patternLists.add(patterns);

        List<SiteRateModel> siteRateModels = new ArrayList<SiteRateModel>();
        siteRateModels.add(siteRateModel);
        siteRateModels.add(siteRateModel2);

        List<BranchModel> branchModels = new ArrayList<BranchModel>();
        branchModels.add(branchModel);
        branchModels.add(branchModel2);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                branchModels,
                siteRateModels,
                true,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: this is 2x the logLikelihood of the 2nd partition)\n\n");

        System.exit(0);



        //START ADDITIONAL TEST #1 - Guy Baele

        System.out.println("-- Test #1 SiteRateModels -- ");
        //alpha in partition 1 reject followed by alpha in partition 2 reject
        System.out.print("Adjust alpha in partition 1: ");
        siteRateModel.setAlpha(0.4);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("Return alpha in partition 1 to original value: ");
        siteRateModel.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (i.e. reject: OK)\n");

        System.out.print("Adjust alpha in partition 2: ");
        siteRateModel2.setAlpha(0.35);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("Return alpha in partition 2 to original value: ");
        siteRateModel2.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (i.e. reject: OK)\n");

        //alpha in partition 1 accept followed by alpha in partition 2 accept
        System.out.print("Adjust alpha in partition 1: ");
        siteRateModel.setAlpha(0.4);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("Adjust alpha in partition 2: ");
        siteRateModel2.setAlpha(0.35);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: same logLikelihood as only setting alpha in partition 2)");

        System.out.print("Return alpha in partition 1 to original value: ");
        siteRateModel.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: alpha in partition 2 has not been returned to original value yet)");

        System.out.print("Return alpha in partition 2 to original value: ");
        siteRateModel2.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + "\n");

        //adjusting alphas in both partitions without explicitly calling getLogLikelihood() in between
        System.out.print("Adjust both alphas in partitions 1 and 2: ");
        siteRateModel.setAlpha(0.4);
        siteRateModel2.setAlpha(0.35);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("Return alpha in partition 2 to original value: ");
        siteRateModel2.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: alpha in partition 1 has not been returned to original value yet)");

        System.out.print("Return alpha in partition 1 to original value: ");
        siteRateModel.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + "\n\n");

        //END ADDITIONAL TEST - Guy Baele


        //START ADDITIONAL TEST #2 - Guy Baele

        System.out.println("-- Test #2 SiteRateModels -- ");
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        //1 siteRateModel shared across 2 partitions
        siteRateModels = new ArrayList<SiteRateModel>();
        siteRateModels.add(siteRateModel);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                branchModels,
                siteRateModels,
                true,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood + "\n");

        System.out.print("Adjust alpha in shared siteRateModel: ");
        siteRateModel.setAlpha(0.4);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: same logLikelihood as only adjusted alpha for partition 1)");
        siteRateModel.setAlpha(0.5);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + "\n\n");

        //END ADDITIONAL TEST - Guy Baele


        //START ADDITIONAL TEST #3 - Guy Baele

        System.out.println("-- Test #3 SiteRateModels -- ");

        siteRateModel = new GammaSiteRateModel("gammaModel");
        siteRateModel.setSubstitutionModel(hky);
        siteRateModel.setRelativeRateParameter(mu);

        siteRateModel2 = new GammaSiteRateModel("gammaModel2");
        siteRateModel2.setSubstitutionModel(hky2);
        siteRateModel2.setRelativeRateParameter(mu);

        siteRateModels = new ArrayList<SiteRateModel>();
        siteRateModels.add(siteRateModel);
        siteRateModels.add(siteRateModel2);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                branchModels,
                siteRateModels,
                true,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.println("logLikelihood = " + logLikelihood + "\n");

        System.out.print("Adjust kappa in partition 1: ");
        hky.setKappa(5.0);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: logLikelihood has not changed?)");

        System.out.print("Return kappa in partition 1 to original value: ");
        hky.setKappa(1.0);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + "\n");

        System.out.print("Adjust kappa in partition 2: ");
        hky2.setKappa(11.0);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood);

        System.out.print("Return kappa in partition 2 to original value: ");
        hky2.setKappa(10.0);
        logLikelihood = treeDataLikelihood.getLogLikelihood();
        System.out.println("logLikelihood = " + logLikelihood + " (i.e. reject: OK)\n\n");

        //END ADDITIONAL TEST - Guy Baele


        //START ADDITIONAL TEST #4 - Guy Baele

        System.out.println("-- Test #4 SiteRateModels -- ");

        SimpleAlignment secondAlignment = createAlignment(moreSequences, Nucleotides.INSTANCE);
        SitePatterns morePatterns = new SitePatterns(secondAlignment, null, 0, -1, 1, true);

        BeagleDataLikelihoodDelegate dataLikelihoodDelegateOne = new BeagleDataLikelihoodDelegate(
                treeModel,
                patterns,
                branchModel,
                siteRateModel, false,
                PartialsRescalingScheme.NONE,
                false);

        TreeDataLikelihood treeDataLikelihoodOne = new TreeDataLikelihood(
                dataLikelihoodDelegateOne,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihoodOne.getLogLikelihood();

        System.out.println("\nBeagleDataLikelihoodDelegate logLikelihood partition 1 (kappa = 1) = " + logLikelihood);

        hky.setKappa(10.0);

        logLikelihood = treeDataLikelihoodOne.getLogLikelihood();

        System.out.println("BeagleDataLikelihoodDelegate logLikelihood partition 1 (kappa = 10) = " + logLikelihood);

        hky.setKappa(1.0);

        BeagleDataLikelihoodDelegate dataLikelihoodDelegateTwo = new BeagleDataLikelihoodDelegate(
                treeModel,
                morePatterns,
                branchModel2,
                siteRateModel2, false,
                PartialsRescalingScheme.NONE,
                false);

        TreeDataLikelihood treeDataLikelihoodTwo = new TreeDataLikelihood(
                dataLikelihoodDelegateTwo,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihoodTwo.getLogLikelihood();

        System.out.println("BeagleDataLikelihoodDelegate logLikelihood partition 2 (kappa = 10) = " + logLikelihood + "\n");

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)patterns),
                Collections.singletonList((BranchModel)branchModel),
                Collections.singletonList((SiteRateModel)siteRateModel),
                true,
                PartialsRescalingScheme.NONE,
                false);


        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.print("Test MultiPartitionDataLikelihoodDelegate 1st partition (kappa = 1):");
        System.out.println("logLikelihood = " + logLikelihood);

        hky.setKappa(10.0);
        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.print("Test MultiPartitionDataLikelihoodDelegate 1st partition (kappa = 10):");
        System.out.println("logLikelihood = " + logLikelihood);

        hky.setKappa(1.0);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                Collections.singletonList((PatternList)morePatterns),
                Collections.singletonList((BranchModel)branchModel2),
                Collections.singletonList((SiteRateModel)siteRateModel2),
                true,
                PartialsRescalingScheme.NONE,
                false);


        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.print("Test MultiPartitionDataLikelihoodDelegate 2nd partition (kappa = 10):");
        System.out.println("logLikelihood = " + logLikelihood + "\n");

        patternLists = new ArrayList<PatternList>();
        patternLists.add(patterns);
        patternLists.add(morePatterns);

        multiPartitionDataLikelihoodDelegate = new MultiPartitionDataLikelihoodDelegate(
                treeModel,
                patternLists,
                branchModels,
                siteRateModels,
                true,
                PartialsRescalingScheme.NONE,
                false);

        treeDataLikelihood = new TreeDataLikelihood(
                multiPartitionDataLikelihoodDelegate,
                treeModel,
                branchRateModel);

        logLikelihood = treeDataLikelihood.getLogLikelihood();

        System.out.print("Test MultiPartitionDataLikelihoodDelegate 2 partitions (kappa = 1, 10): ");

        System.out.println("logLikelihood = " + logLikelihood + " (NOT OK: should be the sum of both separate logLikelihoods)\nKappa value of partition 2 is used to compute logLikelihood for both partitions?");

        //END ADDITIONAL TEST - Guy Baele

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
