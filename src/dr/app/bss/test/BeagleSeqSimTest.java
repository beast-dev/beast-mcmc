/*
 * BeagleSeqSimTest.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.bss.test;

import dr.app.beagle.evomodel.branchmodel.HomogeneousBranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.HKY;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.app.bss.Utils;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ExponentialGrowth;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.ArrayList;

public class BeagleSeqSimTest {

    public static final boolean simulateInPar = true;

    public static void main(String[] args) {

    	
		// start timing
		long tic = System.currentTimeMillis();
        int N = 100; //0000;
		for (int i = 0; i < N; i++) {
//			simulateTopology();
//			simulateOnePartition();
//			simulateTwoPartitions();
			simulateThreePartitions(i, N);
		}
		long toc = System.currentTimeMillis();
		
		long time = toc-tic;
		System.out.println("Time: " + time);
				
		
    } // END: main

    static void simulateTopology() {

        try {

            System.out.println("Test case 1: simulateTopology");

            MathUtils.setSeed(666);

            int sequenceLength = 10;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            Tree tree = importer.importTree(null);

            // set demographic function
            ExponentialGrowth exponentialGrowth = new ExponentialGrowth(
                    Units.Type.YEARS);
            exponentialGrowth.setN0(10);
            exponentialGrowth.setGrowthRate(0.5);

            Taxa taxa = new Taxa();
            for (Taxon taxon : tree.asList()) {

                double absoluteHeight = Utils.getAbsoluteTaxonHeight(taxon,
                        tree);

                taxon.setAttribute(Utils.ABSOLUTE_HEIGHT, absoluteHeight);

                // taxon.setAttribute("date", new Date(absoluteHeight,
                // Units.Type.YEARS, true));

                taxa.addTaxon(taxon);

            }// END: taxon loop

            CoalescentSimulator topologySimulator = new CoalescentSimulator();

            TreeModel treeModel = new TreeModel(topologySimulator.simulateTree(
                    taxa, exponentialGrowth));

            System.out.println(treeModel.toString());

            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25,
                    0.25, 0.25});
            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                    freqs);

            // create substitution model
            Parameter kappa = new Parameter.Default(1, 10);
            HKY hky = new HKY(kappa, freqModel);
            HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
                    hky);

            // create site model
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            // create branch rate model
            BranchRateModel branchRateModel = new DefaultBranchRateModel();

            // create partition
            Partition partition1 = new Partition(treeModel, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    sequenceLength - 1, // to
                    1 // every
            );

            partitionsList.add(partition1);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
                    partitionsList
            );

            System.out.println(simulator.simulate(simulateInPar).toString());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block

    }// END: simulate topology

    static void simulateOnePartition() {

        try {

            MathUtils.setSeed(666);

            System.out.println("Test case 2: simulateOnePartition");

            int sequenceLength = 10;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            Tree tree = importer.importTree(null);
            TreeModel treeModel = new TreeModel(tree);

            // create Frequency Model
            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25,
                    0.25, 0.25});
            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                    freqs);

            // create substitution model
            Parameter kappa = new Parameter.Default(1, 10);
            HKY hky = new HKY(kappa, freqModel);
            HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
                    hky);

            // create site model
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            // create branch rate model
            BranchRateModel branchRateModel = new DefaultBranchRateModel();

            // create partition
            Partition partition1 = new Partition(treeModel, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    sequenceLength - 1, // to
                    1 // every
            );

			Sequence ancestralSequence = new Sequence();
			ancestralSequence.appendSequenceString("TCAAGTGAGG");
			partition1.setAncestralSequence(ancestralSequence);

            partitionsList.add(partition1);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
                    partitionsList
                    // , sequenceLength
            );

            System.out.println(simulator.simulate(simulateInPar).toString());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block

    }// END: simulateOnePartition

    static void simulateTwoPartitions() {

        try {

            System.out.println("Test case 3: simulateTwoPartitions");

            MathUtils.setSeed(666);
            
            int sequenceLength = 11;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            Tree tree = importer.importTree(null);
            TreeModel treeModel = new TreeModel(tree);

            // create Frequency Model
            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25,
                    0.25, 0.25});
            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                    freqs);

            // create substitution model
            Parameter kappa = new Parameter.Default(1, 10);
            HKY hky = new HKY(kappa, freqModel);
            HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
                    hky);

            // create site model
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            // create branch rate model
            BranchRateModel branchRateModel = new DefaultBranchRateModel();

            // create partition
            Partition partition1 = new Partition(treeModel, //
                    substitutionModel, //
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    3, // to
                    1 // every
            );

            // create partition
            Partition Partition = new Partition(treeModel, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    4, // from
                    sequenceLength - 1, // to
                    1 // every
            );

             Sequence ancestralSequence = new Sequence();
             ancestralSequence.appendSequenceString("TCAAGTG");
             Partition.setAncestralSequence(ancestralSequence);

            partitionsList.add(partition1);
            partitionsList.add(Partition);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
                    partitionsList);

            System.out.println(simulator.simulate(simulateInPar).toString());

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block

    }// END: simulateTwoPartitions

    static void simulateThreePartitions(int i, int N) {

        try {
/*
        	>SimSeq1   
        	GCAGGTCTCT
        	>SimSeq2   
        	CATCGGAACT
        	>SimSeq3   
        	CACCTTTACT
*/
            MathUtils.setSeed(666);

//            System.out.println("Test case 4: simulateThreePartitions");

            int sequenceLength = 100000;
            ArrayList<Partition> partitionsList = new ArrayList<Partition>();

            // create tree
            NewickImporter importer = new NewickImporter(
                    "(SimSeq1:73.7468,(SimSeq2:25.256989999999995,SimSeq3:45.256989999999995):18.48981);");
            Tree tree = importer.importTree(null);
            TreeModel treeModel = new TreeModel(tree);

            // create Frequency Model
            Parameter freqs = new Parameter.Default(new double[]{0.25, 0.25,
                    0.25, 0.25});
            FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,
                    freqs);

            // create substitution model
            Parameter kappa = new Parameter.Default(1, 10);
            HKY hky = new HKY(kappa, freqModel);
            HomogeneousBranchModel substitutionModel = new HomogeneousBranchModel(
                    hky);

            // create site model
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel(
                    "siteModel");

            // create branch rate model
            BranchRateModel branchRateModel = new DefaultBranchRateModel();

            // create partition
            Partition partition1 = new Partition(treeModel, //
                    substitutionModel, //
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    0, // from
                    sequenceLength - 1, // to
                    3 // every
            );

            // create partition
            Partition Partition = new Partition(treeModel, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    1, // from
                    sequenceLength - 1, // to
                    3 // every
            );

            // create partition
            Partition partition3 = new Partition(treeModel, //
                    substitutionModel,//
                    siteRateModel, //
                    branchRateModel, //
                    freqModel, //
                    2, // from
                    sequenceLength - 1, // to
                    3 // every
            );

            partitionsList.add(partition1);
            partitionsList.add(Partition);
            partitionsList.add(partition3);

            // feed to sequence simulator and generate data
            BeagleSequenceSimulator simulator = new BeagleSequenceSimulator(
                    partitionsList);

			if (i == (N - 1)) {
				System.out.println(simulator.simulate(simulateInPar).toString());
			} else {
				simulator.simulate(simulateInPar);
			}
            	
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        } // END: try-catch block

    }// END: simulateThreePartitions

}// END: class
