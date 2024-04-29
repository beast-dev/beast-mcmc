/*
 * BASTADensityTester.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SimpleSiteList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.*;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.ScaleOperator;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Guy Baele
 */
public class BASTATreeMoveTester {

    private static final String DEME = "deme";

    private static TreeModel createSpecifiedTree(String t) throws Exception {
        NewickImporter importer = new NewickImporter(t);
        Tree tree = importer.importTree(null);

        return new DefaultTreeModel(tree);//treeModel
    }

    static private String sequences[][] = {
            {"AB192965Japan2004", "AF071228Spain1997", "AF105975Portugal1995", "AF189018Indonesia2005"},
            {
                    "ATGTGGGATCCACTTCTAAATGAATTTCCTGAATCTGTTCACGGATTTCGTTGTATGTTAGCTATTAAATATTTGCAGTCCGTTGAGGAAACTTACGAGCCCAATACATTGGGCCACGATTTAATTAGGGATCTTATATCTGTTGTAAGGGCCCGTGACTATGTCGAAGCGACCAGGCGATATAATCATTTCCACGCCCGCCTCGAAGGTTCGCCGAAGGCTGAACTTCGACAGCCCATACAGCAGCCGTGCTGCTGTCCCCATTGTCCAAGGCACAA---ACAAGCGACGATCATGGACGTACAGGCCCATGTACCGGAAGCCCAGAATATACAGAATGTATCGAAGCCCTGATGTTCCTCGTGGATGTGAAGGCCCATGTAAAGTCCAGTCTTATGAGCAACGGGATGATATTAAGCATACTGGTATTGTTCGTTGTGTTAGTGATGTTACTCGTGGATCTGGAATTACTCACAGAGTGGGTAAGAGGTTCTGTGTTAAATCGATATATTTTTTAGGTAAAGTCTGGATGGATGAAAATATCAAGAAGCAGAATCACACTAATCAGGTCATGTTCTTCTTGGTTCGTGATAGAAGGCCCTATGGAAACAGCCCAATGGATTTTGGACAGGTTTTTAATATGTTCGATAATGAGCCCAGTACCGCAACTGTGAAGAATGATTTGCGGGATAGGTTTCAAGTGATGAGGAAATTTCATGCTACAGTTATTGGTGGGCCCTCTGGAATGAAGGAACAGGCATTAGTTAAGAGATTTTTTAGAATTAACAGTCATGTAACTT-----ATAATCAT-CAGGAGGCAGCCAAGTATGAGAACCATACTGAAAACGCCTTGTTATTGTATATGGCATGTACGCATGCCTCTAATCCAGTGTATGCAACTATGAAAATACGCATCTATTTCTATGATTCAATATCAAATTAA",
                    "ATGTGGGATCCACTTCTAAATGAATTTCCTGAATCTGTTCACGGATTTCGTTGTATGTTAGCTATTAAATATTTGCAGTCCGTTGAGGAAACTTACGAGCCCAATACATTGGGCCACGATTTAATTAGGGATCTTATATCTGTTGTAAGGGCCCGTGACTATGTCGAAGCGACCAGGCGATATAATCATTTCCACGCCCGCCTCGAAGGTTCGCCGAAGGCTGAACTTCGACAGCCCATACAGCAGCCGTGCTGCTGTCCCCATTGTCCAAGGCACAA---ACAAGCGACGATCATGGACGTACAGGCCCATGTACCGGAAGCCCAGAATATACAGAATGTATCGAAGCACTGATGTTCCCCGTGGATGTGAAGGCCCATGTAAAGTACAGTCTTATGAGCAACGGGATGATATTAAGCATACTGGTATTGTTCGTTGTGTTAGTGATGTTACTCGTGGATCTGGAATTACTCACAGAGTGGGTAAGAGGTTCTGTGTTAAATCGATATATTTTTTAGGTAAAGTCTGGATGGATGAAAATATCAAGAAGCAGAACCACACTAATCAGGTCATGTTCTTCTTGGTCCGTGATAGAAGGCCCTATGGAAACAGCCCAATGGATTTTGGACAGGTTTTTAATATGTTCGATAATGAGCCCAGTACCGCAACAGTGAAGAATGATTTGCGGGATAGGTTTCAAGTGATGAGGAAATTTCATGCTACAGTTATTGGTGGGCCCTCTGGAATGAAGGAACAGGCATTAGTTAAGAGGTTTTTTAGAATTAACAGTCATGTAACTT-----ATAATCAT-CAGGAGGCAGCCAAGTACGAGAACCATACTGAAAACGCCTTGTTATTGTATATGGCATGTACGCATGCCTCTAATCCAGTGTATGCTACTATGAAAATACGCATCTATTTCTATGATTCAATATCAAATTAA",
                    "ATGTGGGATCCACTTCTAAATGAATTTCCTGAATCTGTTCACGGATTTCGTTGTATGTTAGCTATTAAATATTTGCAGTCCGTTGAGGAAACTTACGAGCCCAATACATTGGGCCACGATTTAATTAGGGATCTTATATCTGTTGTAAGGGCCCGTGACTATGTCGAAGCGACCAGGCGATATAATCATTTCCACGCCCGCCTCGAAGGTTCGCCGAAGGCTGAACTTCGACAGCCCATACAGCAGCCGTGCTGCTGTCCCCATTGTCCAAGGCACAA---ACAAGCGACGATCATGGACGTACAGGCCCATGTACCGGAAGCCCAGAATATACAGAATGTATCGAAGCCCTGATGTTCCCCGTGGATGTGAAGGCCCATGTAAAGTACAGTCTTATGAGCAACGGGATGATATTAAGCATACTGGTATTGTTCGTTGTGTTAGTGATGTTACTCGTGGATCTGGAATTACTCACAGAGTGGGTAAGAGGTTCTGTGTTAAATCGATATATTTTTTAGGTAAAGTCTGGATGGATGAAAATATCAAGAAGCAGAATCACACTAATCAGGTCATGTTCTTCTTGGTCCGTGATAGAAGGCCCTATGGAAACAGCCCAATGGATTTTGGACAGGTTTTTAATATGTTCGATAATGAGCCCAGTACCGCAACAGTGAAGAATGATTTGCGGGATAGGTTTCAAGTGATGAGGAAATTTCATGCTACAGTTATTGGTGGGCCCTCTGGAATGAAGGAACAGGCATTAGTTAAGAGATTTTTTAGAATTAACAGTCATGTAACTT-----ATAATCAT-CAGGAGGCAGCCAAGTACGAGAACCATACTGAAAACGCCTTGTTATTGTATATGGCATGTACGCATGCCTCTAATCCAGTGTATGCAACTATGAAAATACGCATCTATTTCTATGATTCAATATCAAATTAA",
                    "ATGTGGGATCCTCTTTTAAACGAATTTCCGGATTCTGTTCACGGTTTTCGGTGTATGCTCGCAATAAAGTATTTGCAAGGCGTTGAAGCAACCTACGCCCCTGATACTGTCGGTTACGACCTAGTTCGAGATCTGATCTCAGTTGTTCGTGCGAGCAATTATGCTGAAGCGTGCCGGAGATATAGCCTTTTCCGGTCCCGTATCGAAAGTACGCCGTCGTCTCAATTACGACAGCCCAGGTACCAGCCGTGCTGCTGTACTCACTGCCCTCGGCATAA---ATCGAAAGAAGTCTTGGACTTCTCGGCCTATGTACCGGAAGCCCAGGATTTACCGGATGTACCGAACAGCTGATGTCCCTAGGGGATGTGAAGGTCCTTGCAAGATTCAATCCTTTGAATCTCGACATGATATTGCTCACACCGGTAAGGTTATGTGTGTGACGGATGTTACTCGTGGCGGTGGTTTAACCCACCGTACTGGGAAGAGATTTTGCGTTAAGTCCCTCTATATCCTTGGCAAAATCTGGATGGATGAAAATATCAAGACTAAGAATCACACTAACACGGTCATGTTCTATGTTGTTCGGGATCGTAGACCCTATGGTA-C-TCCT-CAAGATTTTGGACAAGTGTTTAACATGTTCGATAACGAACCTAGCACTGCAACTGTCAAGAATGATCTTCGAGATCGGTTTCAAGTTTTGCGGAAGTTCACGGCAACTGTTGTTGGTGGTCAGTATGCTTGTAAGGAACAAACGTTAGTTAGGAAGTTCATGAGATTGAACAATTATGTTGTTT-----ACAACCAT-CAGGAAACCGCAAAATATGAGAATCATACAGAGAATGCTCTGTTATTGTACATGGCATCTACGCATGCCTCTAACCCTGTGTATGCAACTTTGAAAATTCGGATCTATTTTTATGATTCGATAACAAATTAA"
            }
    };

    private static SimpleAlignment createAlignment(Object[][] taxa_sequence, DataType dataType) {

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);

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

    public static void main(String[] args) {

        MathUtils.setSeed(666);

        // turn off logging to avoid screen noise...
        Logger logger = Logger.getLogger("dr");
        logger.setUseParentHandlers(false);

        System.out.println("EXAMPLE 1: 4 taxa with 2 demes");

        SimpleAlignment alignment = createAlignment(sequences, Nucleotides.INSTANCE);

        TreeModel treeModel;

        try {
            treeModel = createSpecifiedTree("(((AB192965Japan2004:0.42343888910376215,AF189018Indonesia2005:1.4234388891037622)5:9.725099517053918,AF071228Spain1997:3.14853840615768)6:1.7275747971782511,AF105975Portugal1995:2.8761132033359313)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        List<String> states = new ArrayList<String>();
        states.add("Asia");
        states.add("West_Medit");

        GeneralDataType dataType = new GeneralDataType(states);

        Taxa taxonList = new Taxa();
        Taxon taxonOne = new Taxon(sequences[0][0]);
        taxonOne.setAttribute(DEME, states.get(0));
        taxonList.addTaxon(taxonOne);
        Taxon taxonTwo = new Taxon(sequences[0][1]);
        taxonTwo.setAttribute(DEME, states.get(1));
        taxonList.addTaxon(taxonTwo);
        Taxon taxonThree = new Taxon(sequences[0][2]);
        taxonThree.setAttribute(DEME, states.get(1));
        taxonList.addTaxon(taxonThree);
        Taxon taxonFour = new Taxon(sequences[0][3]);
        taxonFour.setAttribute(DEME, states.get(0));
        taxonList.addTaxon(taxonFour);

        int[] pattern = new int[sequences[0].length];
        pattern[0] = 0;
        pattern[1] = 1;
        pattern[2] = 1;
        pattern[3] = 0;

        SimpleSiteList patterns = new SimpleSiteList(dataType, taxonList);
        patterns.addPattern(pattern);

        double[] freqValues = new double[2];
        freqValues[0] = 0.5;
        freqValues[1] = 0.5;
        FrequencyModel freqModel = new FrequencyModel(dataType, freqValues);

        Parameter ratesParameter = new Parameter.Default("migration rates", 2);
        ratesParameter.setParameterValue(0, 1.0);
        ratesParameter.setParameterValue(1, 1.0);

        Parameter popSizesParameter = new Parameter.Default("population sizes", 2);
        popSizesParameter.setParameterValue(0, 1.0);
        popSizesParameter.setParameterValue(1, 1.0);

        //GeneralSubstitutionModel migrationModel = new GeneralSubstitutionModel("migrationModel", dataType, null, ratesParameter, 0);
        SVSComplexSubstitutionModel migrationModel = new SVSComplexSubstitutionModel("migrationModel", dataType, freqModel, ratesParameter, null);

        TaxonList includeSubtree = null;
        List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();

        int subIntervals = 2;

        BranchRateModel branchRateModel = new DefaultBranchRateModel();

        StructuredCoalescentLikelihood structured = null;

        try {
            structured = new StructuredCoalescentLikelihood(treeModel, branchRateModel, popSizesParameter, patterns,null, "", migrationModel, subIntervals, includeSubtree, excludeSubtrees, false);
        } catch (TreeUtils.MissingTaxonException missing) {
            System.out.println("Error thrown in test class dr.evomodel.coalescent.basta.SCLikelihoodTester: " + missing);
        }

        System.out.println("Structured coalescent lnL = " + structured.getLogLikelihood() + "\n");


        //perform scale operator on the root height
        ScaleOperator rootHeightOperator = new ScaleOperator(treeModel.getVariable(0), 0.75);
        rootHeightOperator.doOperation();

        System.out.println("Structured coalescent lnL = " + structured.getLogLikelihood() + "\n");


        //calculate structured coalescent density from scratch to compare result of tree move against
        System.out.println("MODIFIED EXAMPLE 1: 4 taxa with 2 demes (altered root height)");

        try {
            treeModel = createSpecifiedTree("(((AB192965Japan2004:0.42343888910376215,AF189018Indonesia2005:1.4234388891037622):9.725099517053918,AF071228Spain1997:3.14853840615768):1.9314982732284278,AF105975Portugal1995:3.080036679386108)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        try {
            structured = new StructuredCoalescentLikelihood(treeModel, branchRateModel, popSizesParameter, patterns,null, "", migrationModel, subIntervals, includeSubtree, excludeSubtrees, false);
        } catch (TreeUtils.MissingTaxonException missing) {
            System.out.println("Error thrown in test class dr.evomodel.coalescent.basta.SCLikelihoodTester: " + missing);
        }

        System.out.println("Check structured coalescent lnL = " + structured.getLogLikelihood() + "\n");



        //perform scale operator on a different node height
        ScaleOperator nodeHeightOperator = new ScaleOperator(treeModel.getVariable(1), 0.75);
        nodeHeightOperator.doOperation();

        System.out.println("Structured coalescent lnL = " + structured.getLogLikelihood() + "\n");


        //calculate structured coalescent density from scratch to compare result of tree move against
        System.out.println("MODIFIED EXAMPLE 1: 4 taxa with 2 demes (altered node height)");

        try {
            treeModel = createSpecifiedTree("(((AB192965Japan2004:0.42343888910376215,AF189018Indonesia2005:1.4234388891037622):8.95084315354517,AF071228Spain1997:2.3742820426489324):2.7057546367371756,AF105975Portugal1995:3.080036679386108)");
        } catch (Exception e) {
            throw new RuntimeException("Unable to parse Newick tree");
        }

        try {
            structured = new StructuredCoalescentLikelihood(treeModel, branchRateModel, popSizesParameter, patterns,null, "", migrationModel, subIntervals, includeSubtree, excludeSubtrees, false);
        } catch (TreeUtils.MissingTaxonException missing) {
            System.out.println("Error thrown in test class dr.evomodel.coalescent.basta.SCLikelihoodTester: " + missing);
        }

        System.out.println("Check structured coalescent lnL = " + structured.getLogLikelihood() + "\n");

    }

}
