/*
 * AncestralSequenceAnnotator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.tools;


import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodelxml.siteratemodel.GammaSiteModelParser;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.aminoacid.JTT;
import dr.evomodel.substmodel.aminoacid.LG;
import dr.evomodel.substmodel.aminoacid.WAG;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
//import dr.evolution.datatype.AminoAcids;
//import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.*;
import dr.evolution.io.*;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.inference.model.Parameter;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.Version;


import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/*
 * @author Marc A. Suchard
 * @author Wai Lok Sibon Li
 */

public class AncestralSequenceAnnotator {

    private final static Version version = new BeastVersion();

    public final static int MAX_CLADE_CREDIBILITY = 0;
    public final static int MAX_SUM_CLADE_CREDIBILITY = 1;
    public final static int USER_TARGET_TREE = 2;

    public final static int KEEP_HEIGHTS = 0;
    public final static int MEAN_HEIGHTS = 1;
    public final static int MEDIAN_HEIGHTS = 2;

    public final String[] GENERAL_MODELS_LIST = {"EQU"};
    public final String[] NUCLEOTIDE_MODELS_LIST = {"HKY", "TN", "GTR"};
    public final String[] AMINO_ACID_MODELS_LIST = {"JTT1992", "WAG2001", "LG2008", "Empirical\\(.+\\)"};
    public final String[] TRIPLET_MODELS_LIST = {"HKYx3", "TNx3", "GTRx3"};
    public final String[] CODON_MODELS_LIST = {"M0HKY", "M0TN", "M0GTR"};//"M0", "M0\\[.+\\]", ""};

    public AncestralSequenceAnnotator(int burnin,
                                      int heightsOption,
                                      double posteriorLimit,
                                      int targetOption,
                                      String targetTreeFileName,
                                      String inputFileName,
                                      String outputFileName,
                                      String kalignExecutable
    ) throws IOException {

        this.posteriorLimit = posteriorLimit;

        this.kalignExecutable = kalignExecutable;

        attributeNames.add("height");
        attributeNames.add("length");

        System.out.println("Reading trees and simulating internal node states...");

        CladeSystem cladeSystem = new CladeSystem();

        boolean firstTree = true;
        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader);

        TreeExporter exporter;
        if (outputFileName != null)
            exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
        else
            exporter = new NexusExporter(System.out);

        TreeExporter simulationResults = new NexusExporter(new PrintStream(new FileOutputStream(inputFileName + ".out")));
        List<Tree> simulatedTree = new ArrayList<Tree>();

//		burnin = 0;
        //	java.util.logging.Logger.getLogger("dr.evomodel").

        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (firstTree) {
                    Tree unprocessedTree = tree;
                    tree = processTree(tree);
                    setupTreeAttributes(tree);
                    setupAttributes(tree);
                    tree = unprocessedTree;     //This actually does nothing since unprocessedTree was a reference to processedTree in the first place
                    firstTree = false;
                }

                if (totalTrees >= burnin) {
                    addTreeAttributes(tree);
//					Tree savedTree = tree;
                    tree = processTree(tree);
//					System.err.println(Tree.Utils.newick(tree));
                    exporter.exportTree(tree);
//					simulationResults.exportTree(tree);
//					System.exit(-1);
                    simulatedTree.add(tree);
                    cladeSystem.add(tree);

                    totalTreesUsed += 1;
                }
                totalTrees += 1;

            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();

        cladeSystem.calculateCladeCredibilities(totalTreesUsed);

        System.out.println("\tTotal trees read: " + totalTrees);
        if (burnin > 0) {
            System.out.println("\tIgnoring first " + burnin + " trees.");
        }

        simulationResults.exportTrees(simulatedTree.toArray(new Tree[simulatedTree.size()]));

        MutableTree targetTree;

        if (targetOption == USER_TARGET_TREE) {
            if (targetTreeFileName != null) {
                System.out.println("Reading user specified target tree, " + targetTreeFileName);

                importer = new NexusImporter(new FileReader(targetTreeFileName));
                try {
                    targetTree = new FlexibleTree(importer.importNextTree());
                } catch (Importer.ImportException e) {
                    System.err.println("Error Parsing Target Tree: " + e.getMessage());
                    return;
                }
            } else {
                System.err.println("No user target tree specified.");
                return;
            }
        } else if (targetOption == MAX_CLADE_CREDIBILITY) {
            System.out.println("Finding maximum credibility tree...");
            targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName, false));
        } else if (targetOption == MAX_SUM_CLADE_CREDIBILITY) {
            System.out.println("Finding maximum sum clade credibility tree...");
            targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName, true));
        } else {
            throw new RuntimeException("Unknown target tree option");
        }

        System.out.println("Annotating target tree... (this may take a long time)");

//	 	System.out.println("Starting processing....");

//		System.out.println("calling annotateTree");
//		System.out.println("targetTree has "+targetTree.getNodeCount()+" nodes.");
        cladeSystem.annotateTree(targetTree, targetTree.getRoot(), null, heightsOption);

        System.out.println("Processing ended.");

        System.out.println("Writing annotated tree....");
        if (outputFileName != null) {
            exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
            exporter.exportTree(targetTree);
        } else {
            exporter = new NexusExporter(System.out);
            exporter.exportTree(targetTree);
        }


    }



    private abstract class SubstitutionModelLoader {


        public final char[] AA_ORDER = AminoAcids.AMINOACID_CHARS;
        public final char[] NUCLEOTIDE_ORDER = Nucleotides.NUCLEOTIDE_CHARS;
        public final String[] CODON_ORDER = Codons.CODON_TRIPLETS;

        protected SubstitutionModel substModel;
        protected FrequencyModel freqModel;
        private String modelType;
        protected String substModelName;
        protected String[] charList;
        protected DataType dataType;

        SubstitutionModelLoader(Tree tree, String modelType, DataType dataType) {
            this.dataType = dataType;
            this.modelType = modelType;
            load(tree, modelType);
        }

        /* An artifact of old code? */
        SubstitutionModelLoader(String name) {
            this.substModelName = name;
        }

        public SubstitutionModel getSubstitutionModel() {
            return substModel;
        }

        public FrequencyModel getFrequencyModel() {
            return freqModel;
        }

        public String getModelType() {
            return modelType;
        }

        public String getSubstModel() {
            return substModelName;
        }

        public String[] getCharList() {
            return charList;
        }
        public void setCharList(String[] cl) {
            System.arraycopy(cl, 0, charList, 0, cl.length);
            //charList = cl.clone();
        }

        public abstract DataType getDataType();

        protected abstract void modelSpecifics(Tree tree, String modelType);

        public void load(Tree tree, String modelType) {

            substModelName = modelType.replaceFirst("\\*.+","").replaceFirst("\\+.+","").trim();

            loadFrequencyModel(tree);

            modelSpecifics(tree, modelType);
            printLogLikelihood(tree);
        }


        private void loadFrequencyModel(Tree tree) {
            final String[] AA_ORDER = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R",
                    "S","T","V","W","Y"}; //"ACDEFGHIKLMNPQRSTVWY".split("");//AminoAcids.AMINOACID_CHARS;
            final String[] DNA_NUCLEOTIDE_ORDER = {"A", "C", "G", "T"};//"ACGT".split(""); //Nucleotides.NUCLEOTIDE_CHARS;
            final String[] RNA_NUCLEOTIDE_ORDER = {"A", "C", "G", "U"};//"ACGU".split(""); //Nucleotides.NUCLEOTIDE_CHARS;
//            final String[] CODON_ORDER = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
//                    "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
//                    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
//                    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
//                    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
//                    "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
//                    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
//                    "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};// Codons.CODON_TRIPLETS;
            final String[] DNA_CODON_ORDER = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
                "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
                "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                /*"TAA",*/ "TAC", /*"TAG",*/ "TAT", "TCA", "TCC", "TCG", "TCT", /* Minus the stop and start codons */
                /*"TGA",*/ "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};// Codons.CODON_TRIPLETS;

            final String[] RNA_CODON_ORDER = {"AAA", "AAC", "AAG", "AAU", "ACA", "ACC", "ACG", "ACU",
                "AGA", "AGC", "AGG", "AGU", "AUA", "AUC", "AUG", "AUU",
                "CAA", "CAC", "CAG", "CAU", "CCA", "CCC", "CCG", "CCU",
                "CGA", "CGC", "CGG", "CGU", "CUA", "CUC", "CUG", "CUU",
                "GAA", "GAC", "GAG", "GAU", "GCA", "GCC", "GCG", "GCU",
                "GGA", "GGC", "GGG", "GGU", "GUA", "GUC", "GUG", "GUU",
                /*"UAA",*/ "UAC", /*"UAG",*/ "UAU", "UCA", "UCC", "UCG", "UCU", /* Minus the stop and start codons */
                /*"UGA",*/ "UGC", "UGG", "UGU", "UUA", "UUC", "UUG", "UUU"};// Codons.CODON_TRIPLETS;

            /* For BAli-Phy, even if F=constant, you can still extract the frequencies this way. */

            /* Obtain the equilibrium base frequencies for the model */
            double[] freq = new double[0];
            String[] charOrder = new String[0];

            if(getDataType().getClass().equals(GeneralDataType.class)) {
                ArrayList<String> tempCharOrder = new ArrayList<String>(freq.length);
                for (Iterator<String> i = tree.getAttributeNames(); i.hasNext();) {
                    String name = i.next();
                    if (name.startsWith("pi")) { /* the pi in the output files contains the frequencies */
                        String character = name.substring(2, name.length());
                        tempCharOrder.add(character);
                    }
                }
                charOrder = tempCharOrder.toArray(new String[tempCharOrder.size()]);
                //this.charList = tempCharOrder.toArray(new String[tempCharOrder.size()]);
                freq = new double[charOrder.length];

            }
            else if(getDataType().getClass().equals(Nucleotides.class)) {
                if(tree.getAttribute("piT") != null) {
                    charOrder = DNA_NUCLEOTIDE_ORDER;
                    freq = new double[charOrder.length];
                }
                else if(tree.getAttribute("piU") != null) {
                    charOrder = RNA_NUCLEOTIDE_ORDER;
                    freq = new double[charOrder.length];
                }
                else {
                    throw new RuntimeException("Not proper nucleotide data");
                }
            }
            else if(getDataType().getClass().equals(AminoAcids.class)) {

                charOrder = AA_ORDER;
                freq = new double[charOrder.length];
            }

            else if(getDataType().getClass().equals(Codons.class)) {
                if(tree.getAttribute("piAAT") != null) {
                    charOrder = DNA_CODON_ORDER;
                    freq = new double[charOrder.length];
                }
                else if(tree.getAttribute("piAAU") != null) {
                    charOrder = RNA_CODON_ORDER;
                    freq = new double[charOrder.length];
                }
                else{
                    throw new RuntimeException("Base frequencies do not fit those for a codon model or not proper nucleotide data for codons\n" +
                        "If you are using F=nucleotides models for codon models, they are currently not supported in BEAST");
                }

            }
            else {
                throw new RuntimeException("Datatype unknown! (This error message should never be seen, contact Sibon)");
            }
    //            // This if statement is reserved for triplet data
    //            else if(getDataType().getClass().equals(Nucleotides.class)) {
    //                // This is wrong because for triplets
    //                freq = new double[Nucleotides.NUCLEOTIDE_CHARS.length];
    //
    //            }

            int cnt = 0;
            double sum = 0;
            //charList = "";
            ArrayList<String> tempCharList = new ArrayList<String>(freq.length);
            for (Iterator<String> i = tree.getAttributeNames(); i.hasNext();) {
                String name = i.next();
                if (name.startsWith("pi")) { /* the pi in the output files contains the frequencies */
                    String character = name.substring(2, name.length());
                    tempCharList.add(character);
                    //charList = charList.concat(character);
                    Double value = (Double) tree.getAttribute(name);
                    freq[cnt++] = value;
                    sum += value;
                }
            }
            charList = tempCharList.toArray(new String[tempCharList.size()]);

//            for(int j=0; j<charList.length; j++) {
//                System.out.println("charizard lists " + charList[j]);
//            }

            /* Order the frequencies correctly */
            double[] freqOrdered = new double[freq.length];
            for (int i = 0; i < freqOrdered.length; i++) {
                int index = -5;
                search: for (int j = 0; j < charList.length; j++) {
                    if(charList[j].equals(charOrder[i])) {
                        index = j;
                        break search;
                    }
                }

                freqOrdered[i] = freq[index] / sum;
                //System.out.println(" no fried " + freqOrdered.length + "\t" + freq.length + "\t" + charOrder[i] + "\t" + freqOrdered[i]); //ddd
            }
            this.freqModel = new FrequencyModel(getDataType(), new Parameter.Default(freqOrdered));

        }

        protected boolean doPrint = false;
        private void printLogLikelihood(Tree tree) {
            if (doPrint) {
                Double logLikelihood = Double.parseDouble(tree.getAttribute(LIKELIHOOD).toString());
                if (logLikelihood != null)
                    System.err.printf("%5.1f", logLikelihood);
            }
        }
    }

    private class GeneralSubstitutionModelLoader extends SubstitutionModelLoader {

        private final String EQU_TEXT = "EQU";


        private GeneralSubstitutionModelLoader(Tree tree, String modelType) {
            super(tree, modelType, new GeneralDataType(new String[0]));
            setGeneralDataType();
            throw new RuntimeException("General substitution model is currently not stable and should not be used");
        }
        protected void modelSpecifics(Tree tree, String modelType) {
            if(substModelName.equals(EQU_TEXT)) {
                if(freqModel.getFrequencyCount() != charList.length) {
                    System.err.println("Frequency model length does not match character list length, " +
                            "GeneralSubstitutionModelLoader");
                    System.exit(-1);
                }

                /* Equivalent to a JC model but for all states */
                //TODO CHECK IF THIS IS CORRECT
                double[] rates = new double[(charList.length * (charList.length - 1)) / 2];
                for(int i=0; i<rates.length; i++) {
                    rates[i] = 1.0;
                }
                System.out.println("Number of site transition rate categories (debuggin): " + rates.length);
                //substModel = new GeneralSubstitutionModel(freqModel.getDataType(), freqModel, new Parameter.Default(rates), 1);
                substModel = new GeneralSubstitutionModel(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, freqModel.getDataType(), freqModel, new Parameter.Default(rates), 1, null);
            }
        }
        public DataType getDataType() {
            setGeneralDataType();
            return dataType;
        }

        public void setGeneralDataType() {
            if(charList!=null && dataType.getStateCount()!=charList.length) {
                GeneralDataType gdt = new GeneralDataType(charList);
                gdt.addAmbiguity("-", charList);
                gdt.addAmbiguity("X", charList);
                gdt.addAmbiguity("?", charList);
                dataType = gdt;
            }
        }
    }

    private class NucleotideSubstitutionModelLoader extends SubstitutionModelLoader {
        protected static final String HKY_TEXT = "HKY";
        protected static final String TN_TEXT = "TN";
        protected static final String GTR_TEXT = "GTR";

        private NucleotideSubstitutionModelLoader(Tree tree, String modelType) {
            super(tree, modelType, Nucleotides.INSTANCE);
        }
        protected void modelSpecifics(Tree tree, String modelType) {
            if(substModelName.equals(HKY_TEXT)) {
                double kappa = Double.parseDouble(tree.getAttribute("HKY_kappa").toString());
                //double kappa = (Double) tree.getAttribute("HKY_kappa");
                //double kappa = (Double) tree.getAttribute("HKY\\:\\:kappa");
                //double kappa = (Double) tree.getAttribute("HKY::kappa");
                //double kappa = (Double) tree.getAttribute("kappa");
                substModel = new HKY(new Parameter.Default(kappa), freqModel);

            }
            if(substModelName.equals(TN_TEXT)) {

                double kappa1 = Double.parseDouble(tree.getAttribute("TN_kappa(pur)").toString());
                double kappa2 = Double.parseDouble(tree.getAttribute("TN_kappa(pyr)").toString());
                //double kappa1 = (Double) tree.getAttribute("TN_kappa(pur)");
                //double kappa2 = (Double) tree.getAttribute("TN_kappa(pyr)");
                System.err.println("Sorry, TN substitution model is not yet implemented in BEAST-Beagle");
                System.exit(0);

                //TODO Tamura-Nei model
                //substModel = new TN93(new Parameter.Default(kappa1), new Parameter.Default(kappa2), freqModel);
            }
            if(substModelName.equals(GTR_TEXT)) {
                /* It should be noted that BAli-Phy uses TC instead of CT and GC instead of CG */
                //double rateACValue = (Double) tree.getAttribute("GTR_AC");
                //double rateAGValue = (Double) tree.getAttribute("GTR_AG");
                //double rateATValue = (Double) tree.getAttribute("GTR_AT");
                //double rateCGValue = (Double) tree.getAttribute("GTR_GC");
                //double rateCTValue = (Double) tree.getAttribute("GTR_TC");
                //double rateGTValue = (Double) tree.getAttribute("GTR_GT");
                double rateACValue = Double.parseDouble(tree.getAttribute("GTR_AC").toString());
                double rateAGValue = Double.parseDouble(tree.getAttribute("GTR_AG").toString());
                double rateATValue = Double.parseDouble(tree.getAttribute("GTR_AT").toString());
                double rateCGValue = Double.parseDouble(tree.getAttribute("GTR_GC").toString());
                double rateCTValue = Double.parseDouble(tree.getAttribute("GTR_TC").toString());
                double rateGTValue = Double.parseDouble(tree.getAttribute("GTR_GT").toString());


                substModel = new GTR(new Parameter.Default(rateACValue), new Parameter.Default(rateAGValue),
                        new Parameter.Default(rateATValue), new Parameter.Default(rateCGValue),
                        new Parameter.Default(rateCTValue), new Parameter.Default(rateGTValue), freqModel);
            }
        }
        public DataType getDataType() {
            return dataType; //Potential
        }
    }

    private class AminoAcidSubstitutionModelLoader extends SubstitutionModelLoader {

        protected final String JTT_TEXT = "JTT1992";
        protected final String WAG_TEXT = "WAG2001";
        protected final String LG_TEXT = "LG2008";
        protected final String Empirical_TEXT = "Empirical(.+).+";

        private AminoAcidSubstitutionModelLoader(Tree tree, String modelType) {
            super(tree, modelType, AminoAcids.INSTANCE);
        }
        protected void modelSpecifics(Tree tree, String modelType) {
            if(substModelName.equals(JTT_TEXT)) {
                substModel = new EmpiricalAminoAcidModel(JTT.INSTANCE, freqModel);
            }
            if(substModelName.equals(WAG_TEXT)) {
                substModel = new EmpiricalAminoAcidModel(WAG.INSTANCE, freqModel);
            }
            if(substModelName.equals(LG_TEXT)) {
                substModel = new EmpiricalAminoAcidModel(LG.INSTANCE, freqModel);
            }
            //todo Allow proper file input of Empirical amino-acid models
            if(substModelName.matches(Empirical_TEXT)) {
                String empiricalModelFileName = substModelName.replaceFirst("Empirical\\(", "").replaceFirst("\\).*","");
                if (empiricalModelFileName.equals("wag.dat")) {
                    substModel = new EmpiricalAminoAcidModel(WAG.INSTANCE, freqModel);
                } else if(empiricalModelFileName.equals("jtt.dat")) {
                    substModel = new EmpiricalAminoAcidModel(JTT.INSTANCE, freqModel);
                } else if(empiricalModelFileName.equals("lg.dat")) {
                    substModel = new EmpiricalAminoAcidModel(LG.INSTANCE, freqModel);
                } else {
                    System.err.println("Sorry, AncestralSequenceAnnotator does not currently support other files");
                    System.err.println("Soon, we will allow users to enter a file");
                    System.exit(0);
                }
            }
        }
        public DataType getDataType() {
            return dataType; //Potential
        }
    }

    private class TripletSubstitutionModelLoader extends SubstitutionModelLoader {

        protected final String HKYx3_TEXT = "HKYx3";
        protected final String TNx3_TEXT = "TNx3";
        protected final String GTRx3_TEXT = "GTRx3";

        private TripletSubstitutionModelLoader(Tree tree, String modelType) {
            super(tree, modelType, Nucleotides.INSTANCE);
        }
        protected void modelSpecifics(Tree tree, String modelType) {
            if(substModelName.equals(HKYx3_TEXT)) {
                System.err.println("Sorry, HKYx3 substitution model is not yet implemented in BEAST");
                System.exit(0);
                //substModel = new HKY();
            }
            if(substModelName.equals(TNx3_TEXT)) {
                System.err.println("Sorry, TNx3 substitution model is not yet implemented in BEAST");
                System.exit(0);
                //substModel = new TN93();
            }
            if(substModelName.equals(GTRx3_TEXT)) {
                System.err.println("Sorry, GTRx3 substitution model is not yet implemented in BEAST");
                System.exit(0);
                //substModel = new GTR();
            }
        }
        public DataType getDataType() {
            return dataType; // Is this right?
        }
    }
    private class CodonSubstitutionModelLoader extends SubstitutionModelLoader {

//        private final String M0_TEXT = "M0"; // Not necessary since this never actually shows up in BAli-Phy output
        protected final String M0_NUC_TEXT = "M0\\w+";
        //private final String HKY_TEXT = "HKY";
        //private final String TN_TEXT = "TN";
        //private final String GTR_TEXT = "GTR";

        private CodonSubstitutionModelLoader(Tree tree, String modelType) {
            super(tree, modelType, Codons.UNIVERSAL); //Potential
        }
        protected void modelSpecifics(Tree tree, String modelType) {
            //String codonNucleotideModel = substModelName.substring(substModelName.indexOf("\\\\[")+1, substModelName.indexOf("\\\\]"));
            String codonNucleotideModel = substModelName.substring(substModelName.indexOf("M0")+2, substModelName.length());


//            if(substModelName.equals(M0_TEXT)) {
//                /* HKY is default */
//                codonNucleotideModel = NucleotideSubstitutionModelLoader.HKY_TEXT;
//
//
//            }
            if(substModelName.matches(M0_NUC_TEXT)) {
                /* M0_omega may be *M0_omega, depending on whether M2 etc. are used */
                double omega = Double.parseDouble(tree.getAttribute("M0_omega").toString());
                //omega = (Double) tree.getAttribute("\\M0_omega");
                if(codonNucleotideModel.equals(NucleotideSubstitutionModelLoader.HKY_TEXT)) {
                    //double kappa = (Double) tree.getAttribute("HKY_kappa");
                    double kappa = Double.parseDouble(tree.getAttribute("HKY_kappa").toString());
                    //substModel = new YangCodonModel(Codons.UNIVERSAL, new Parameter.Default(omega),
                            //new Parameter.Default(kappa), freqModel);
                    substModel = new GY94CodonModel(Codons.UNIVERSAL, new Parameter.Default(omega),
                            new Parameter.Default(kappa), freqModel);
                }
                if(codonNucleotideModel.equals(NucleotideSubstitutionModelLoader.TN_TEXT)) {
                    //double kappa1 = (Double) tree.getAttribute("TN_kappa(pur)");
                    //double kappa2 = (Double) tree.getAttribute("TN_kappa(pyr)");
                    double kappa1 = Double.parseDouble(tree.getAttribute("TN_kappa(pur)").toString());
                    double kappa2 = Double.parseDouble(tree.getAttribute("TN_kappa(pyr)").toString());
                    System.err.println("Sorry, M0[TN] substitution model is not yet implemented in BEAST");
                    System.exit(0);
                }
                if(codonNucleotideModel.equals(NucleotideSubstitutionModelLoader.GTR_TEXT)) {
                    double rateACValue = Double.parseDouble(tree.getAttribute("GTR_AC").toString());
                    double rateAGValue = Double.parseDouble(tree.getAttribute("GTR_AG").toString());
                    double rateATValue = Double.parseDouble(tree.getAttribute("GTR_AT").toString());
                    double rateCGValue = Double.parseDouble(tree.getAttribute("GTR_GC").toString());
                    double rateCTValue = Double.parseDouble(tree.getAttribute("GTR_TC").toString());
                    double rateGTValue = Double.parseDouble(tree.getAttribute("GTR_GT").toString());
                    //double rateACValue = (Double) tree.getAttribute("GTR_AC");
                    //double rateAGValue = (Double) tree.getAttribute("GTR_AG");
                    //double rateATValue = (Double) tree.getAttribute("GTR_AT");
                    //double rateCGValue = (Double) tree.getAttribute("GTR_GC");
                    //double rateCTValue = (Double) tree.getAttribute("GTR_TC");
                    //double rateGTValue = (Double) tree.getAttribute("GTR_GT");
                    System.err.println("Sorry, M0[GTR] substitution model is not yet implemented in BEAST");
                    System.exit(0);
                }
                // If +m2 then *M0
            }
        }
        public DataType getDataType() {
            return dataType; // Is this right too? Just use the universal codon table?
        }
    }

    /*
     * This method is equivalent to the SubstitutionModelLoader without having to
     * be object orientated and can be much more flexible.
     */
    private GammaSiteRateModel loadSiteModel(Tree tree) {

        String modelType = (String) tree.getAttribute(SUBST_MODEL);

        /* Identify the datatype and substitution model. Load the model */
        SubstitutionModelLoader sml = null;
        String substModelName = modelType.replaceFirst("\\*.+","").replaceFirst("\\+.+","").trim();
        System.out.println("Basic Substitution Model is " + substModelName);
        for(int i = 0; i<GENERAL_MODELS_LIST.length; i++) {
            if(substModelName.matches(GENERAL_MODELS_LIST[i])) {
                sml = new GeneralSubstitutionModelLoader(tree, modelType);
            }
        }
        for(int i = 0; i<NUCLEOTIDE_MODELS_LIST.length; i++) {
            if(substModelName.matches(NUCLEOTIDE_MODELS_LIST[i])) {
                sml = new NucleotideSubstitutionModelLoader(tree, modelType);
            }
        }
        for(int i = 0; i<AMINO_ACID_MODELS_LIST.length; i++) {
            if(substModelName.matches(AMINO_ACID_MODELS_LIST[i])) {
                sml = new AminoAcidSubstitutionModelLoader(tree, modelType);
            }
        }
        for(int i = 0; i<TRIPLET_MODELS_LIST.length; i++) {
            if(substModelName.matches(TRIPLET_MODELS_LIST[i])) {
                sml = new TripletSubstitutionModelLoader(tree, modelType);
            }
        }
        for(int i = 0; i<CODON_MODELS_LIST.length; i++) {
            if(substModelName.matches(CODON_MODELS_LIST[i])) {
                sml = new CodonSubstitutionModelLoader(tree, modelType);
            }
        }

        if (sml.getSubstitutionModel() == null) {
            System.err.println("Substitution model type '" + modelType + "' not implemented");
            System.exit(-1);
        }

        //SiteModel siteModel = new GammaSiteModel(sml.getSubstitutionModel(), new Parameter.Default(1.0), null, 0, null);
        //SiteModel siteModel = new GammaSiteModel(sml.getSubstitutionModel(), null, null, 0, null);

        String siteRatesModels = modelType.substring(modelType.indexOf("+")+1, modelType.length());
        //String[] siteRatesModels = siteRatesParameters.split(" + ");
        System.out.println("Site rate models: " + siteRatesModels);
        if(sml.getSubstitutionModel().getDataType().getClass().equals(Codons.class) &&
                siteRatesModels.length() > 0) { /* For codon site models */

            if(siteRatesModels.indexOf("+M2") >= 0) { /* M2 */
                System.out.println("Site model - M2 Codon site model used");
                Parameter m2FrequencyAAInv = new Parameter.Default(Double.parseDouble(tree.getAttribute("M2_f[AA INV]").toString()));
                Parameter m2FrequencyNeutral = new Parameter.Default(Double.parseDouble(tree.getAttribute("M2_f[Neutral]").toString()));
                Parameter m2FrequencySelected = new Parameter.Default(Double.parseDouble(tree.getAttribute("M2_f[Selected]").toString()));
                Parameter m2Omega = new Parameter.Default(Double.parseDouble(tree.getAttribute("M2_omega").toString()));
                //Parameter m2FrequencyAAInv = new Parameter.Default((Double) tree.getAttribute("M2_f[AA INV]"));
                //Parameter m2FrequencyNeutral = new Parameter.Default((Double) tree.getAttribute("M2_f[Neutral]"));
                //Parameter m2FrequencySelected = new Parameter.Default((Double) tree.getAttribute("M2_f[Selected]"));
                //Parameter m2Omega = new Parameter.Default((Double) tree.getAttribute("M2_omega"));
                System.err.println("Sorry, M2 substitution model is not yet implemented in BEAST");
                System.exit(0);
            }
            else if(siteRatesModels.indexOf("+M3") >= 0) { /* M3 */
                System.out.println("Site model - M3 Codon site model used");
                int numberOfBins = Integer.parseInt(siteRatesModels.replaceFirst(".+M3\\[","").replaceFirst("\\].+", ""));
                System.out.println(" + M3 n value: " + numberOfBins);
                Parameter[] m3Frequencies = new Parameter[numberOfBins];
                Parameter[] m3Omegas = new Parameter[numberOfBins];
                for(int i=1; i<=numberOfBins; i++) {
                    m3Frequencies[i-1] = new Parameter.Default(Double.parseDouble(tree.getAttribute("M3_f"+i).toString()));
                    m3Omegas[i-1] = new Parameter.Default(Double.parseDouble(tree.getAttribute("M3_omega"+i).toString()));
                    //m3Frequencies[i-1] = new Parameter.Default((Double) tree.getAttribute("M3_f"+i));
                    //m3Omegas[i-1] = new Parameter.Default((Double) tree.getAttribute("M3_omega"+i));
                }
                System.err.println("Sorry, M3 substitution model is not yet implemented in BEAST");
                System.exit(0);
            }
            else if(siteRatesModels.indexOf("+M0_omega~Beta(") >= 0) { /* M7 */
                System.out.println("Site model - M7 Codon site model used");
                int numberOfBins = Integer.parseInt(siteRatesModels.replaceFirst("M0_omega~Beta\\(","").replaceFirst("\\)", ""));
                System.out.println(" + M7 n value: " + numberOfBins);
                Parameter m7BetaMu = new Parameter.Default(Double.parseDouble(tree.getAttribute("beta_mu").toString()));
                Parameter m7BetaVarMu = new Parameter.Default(Double.parseDouble(tree.getAttribute("beta_Var/mu").toString()));
                //Parameter m7BetaMu = new Parameter.Default((Double) tree.getAttribute("beta_mu"));
                //Parameter m7BetaVarMu = new Parameter.Default((Double) tree.getAttribute("beta_Var/mu"));
                System.err.println("Sorry, M7 substitution model is not yet implemented in BEAST");
                System.exit(0);
            }
        }
        else if(siteRatesModels.length() > 0) { /* i.e. for other data types. */
            /* Do gamma/lognormal + pinv */
            Parameter pInvParameter = null;
            int categories = -1;
            Parameter alphaParameter = null;

            //System.out.println("Greatest story ever told! " + siteRatesModels);

            if(siteRatesModels.indexOf("+INV") >= 0) {
                System.out.println("Site model -  proportion of invariable sites used");
                //pInvParameter = new Parameter.Default(((Double) tree.getAttribute("INV_p")).doubleValue());
                pInvParameter = new Parameter.Default(Double.parseDouble(tree.getAttribute("INV_p").toString()));
            }
            if(siteRatesModels.indexOf("+rate~Gamma(") >= 0) {
                System.out.println("Site model - gamma site rate heterogeneity used");
                categories = Integer.parseInt(siteRatesModels.replaceFirst(".+rate~Gamma\\(", "").replaceFirst("\\).*",""));
                //double sigmaMu = (Double) tree.getAttribute("gamma_sigma/mu");
                double sigmaMu = Double.parseDouble(tree.getAttribute("gamma_sigma/mu").toString());
                sigmaMu = (1.0/sigmaMu) * (1.0/sigmaMu); /* BAli-Phy is parameterised by sigma/mu instead of alpha */
                alphaParameter = new Parameter.Default(sigmaMu);
            }
            else if(siteRatesModels.indexOf("+rate~LogNormal(") >= 0) {
                // TODO implement lognormal site model
                System.out.println("Site model - lognormal site rate heterogeneity used");
                System.err.println("Sorry, lognormal site rates are not yet implemented in BEAST");
                System.exit(0);

                categories = Integer.parseInt(siteRatesModels.replaceFirst(".+rate~LogNormal\\(", "").replaceFirst("\\).*",""));
                //double sigmaMu = (Double) tree.getAttribute("log-normal_sigma/mu");
                double sigmaMu = Double.parseDouble(tree.getAttribute("log-normal_sigma/mu").toString());
                sigmaMu = (1.0/sigmaMu) * (1.0/sigmaMu); /* BAli-Phy is parameterised by sigma/mu instead of alpha */
                alphaParameter = new Parameter.Default(sigmaMu);
            }
            else if(siteRatesModels.indexOf("+GAMMA(") >= 0) { /* For BEAST output */
                System.out.println("Site model - gamma site rate heterogeneity used");
                categories = Integer.parseInt(siteRatesModels.replaceFirst(".+GAMMA\\(", "").replaceFirst("\\).*",""));
                //double sigmaMu = (Double) tree.getAttribute("gamma_sigma/mu");
                double alpha = Double.parseDouble(tree.getAttribute("alpha").toString());
                alphaParameter = new Parameter.Default(alpha);
            }

            //System.out.println("alpha and pinv parameters: " + alphaParameter.getParameterValue(0) + "\t" + pInvParameter.getParameterValue(0));
            //GammaSiteRateModel siteModel = new GammaSiteRateModel(sml.getSubstitutionModel(), new Parameter.Default(1.0), alphaParameter, categories, pInvParameter);
            GammaSiteRateModel siteModel = new GammaSiteRateModel(GammaSiteModelParser.SITE_MODEL, new Parameter.Default(1.0), alphaParameter, categories, pInvParameter);
            siteModel.setSubstitutionModel(sml.getSubstitutionModel());
            //SiteModel siteModel = new GammaSiteModel(sml.getSubstitutionModel(), new Parameter.Default(1.0), new Parameter.Default(1.0), 1, new Parameter.Default(0.5));
            //SiteModel siteModel = new GammaSiteModel(sml.getSubstitutionModel(), null, null, 0, null);
            return siteModel;
        }

        /* Default with no gamma or pinv */
        //SiteRateModel siteModel = new GammaSiteRateModel(sml.getSubstitutionModel());
        GammaSiteRateModel siteModel = new GammaSiteRateModel(GammaSiteModelParser.SITE_MODEL);
        siteModel.setSubstitutionModel(sml.getSubstitutionModel());
        return siteModel;

    }


    public static final String KAPPA_STRING = "kappa";
    //public static final String SEQ_STRING = "states"; // For BEAST input files
    public static String SEQ_STRING = "seq";
    //public static final String SEQ_STRING_2 = "states"; // For BEAST input files
    public static final String NEW_SEQ = "newSeq";
    public static final String TAG = "tag";
    public static final String LIKELIHOOD = "lnL";
    public static final String SUBST_MODEL = "subst";


    //	public static final String WAG_STRING = "Empirical(Data/wag.dat)*pi";
    private final int die = 0;

    private Tree processTree(Tree tree) {

        // Remake tree to fix node ordering - Marc

        GammaSiteRateModel siteModel = loadSiteModel(tree);

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(siteModel.getSubstitutionModel().getDataType());
        if(siteModel.getSubstitutionModel().getDataType().getClass().equals(Codons.class)) {
            //System.out.println("trololo");
            alignment.setDataType(Nucleotides.INSTANCE);
        }
        //System.out.println("BOO BOO " + siteModel.getSubstitutionModel().getDataType().getClass().getName()+"\t" + Codons.UNIVERSAL.getClass().getName() + "\t" + alignment.getDataType().getClass().getName());


        // Get sequences
        String[] sequence = new String[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            sequence[i] = (String) tree.getNodeAttribute(node, SEQ_STRING);
            if (tree.isExternal(node)) {
                Taxon taxon = tree.getNodeTaxon(node);
                alignment.addSequence(new Sequence(taxon, sequence[i]));
                //System.out.println("seq " + sequence[i]);
            }
        }

        // Make evolutionary model

        BranchRateModel rateModel = new StrictClockBranchRates(new Parameter.Default(1.0));
        FlexibleTree flexTree;
        if(siteModel.getSubstitutionModel().getDataType().getClass().equals(Codons.class)) {
            ConvertAlignment convertAlignment  = new ConvertAlignment(siteModel.getSubstitutionModel().getDataType(), ((Codons) siteModel.getSubstitutionModel().getDataType()).getGeneticCode(), alignment);
            flexTree = sampleTree(tree, convertAlignment, siteModel, rateModel);
            //flexTree = sampleTree(tree, alignment, siteModel, rateModel);
        }
        else {
            flexTree = sampleTree(tree, alignment, siteModel, rateModel);
        }
        introduceGaps(flexTree, tree);

        return flexTree;
    }





    private void introduceGaps(FlexibleTree flexTree, Tree gapTree) {
        // I forget what this function was supposed to do. - Marc
    }


    public static final char GAP = '-';

    boolean[] bit = null;

    private FlexibleTree sampleTree(Tree tree, PatternList alignment, GammaSiteRateModel siteModel, BranchRateModel rateModel) {
        FlexibleTree flexTree = new FlexibleTree(tree, true);
        flexTree.adoptTreeModelOrdering();
        FlexibleTree finalTree = new FlexibleTree(tree);
        finalTree.adoptTreeModelOrdering();
        TreeModel treeModel = new TreeModel(tree);

        // Turn off noisy logging by TreeLikelihood constructor
        Logger logger = Logger.getLogger("dr.evomodel");
        boolean useParentHandlers = logger.getUseParentHandlers();
        logger.setUseParentHandlers(false);

//        AncestralStateTreeLikelihood likelihood = new AncestralStateTreeLikelihood(
//                alignment,
//                treeModel,
//                siteModel,
//                rateModel,
//                false, true,
//                alignment.getDataType(),
//                TAG,
//                false);
        AncestralStateBeagleTreeLikelihood likelihood = new AncestralStateBeagleTreeLikelihood(
                alignment,
                treeModel,
                new HomogeneousBranchModel(siteModel.getSubstitutionModel()),
                siteModel,
                rateModel,
                null,
                false,
                PartialsRescalingScheme.DEFAULT,
                true,
                null,
                alignment.getDataType(),
                TAG,
                false, true
                );

//        PatternList patternList, TreeModel treeModel,
//                                              BranchSiteModel branchSiteModel,
// SiteRateModel siteRateModel,
//                                              BranchRateModel branchRateModel,
// boolean useAmbiguities,
//                                              PartialsRescalingScheme scalingScheme,
//                                              Map<Set<String>, Parameter> partialsRestrictions,
//                                              final DataType dataType,
//                                              final String tag,
//                                              SubstitutionModel substModel,
//                                              boolean useMAP,
//                                              boolean returnML) {

//        PatternList patternList, TreeModel treeModel,
//                                        SiteModel siteModel, BranchRateModel branchRateModel,
//                                        boolean useAmbiguities, boolean storePartials,
//                                        final DataType dataType,
//                                        final String tag,
//                                        boolean forceRescaling,
//                                        boolean useMAP,
//                                        boolean returnML) {


        logger.setUseParentHandlers(useParentHandlers);

        // Sample internal nodes
        likelihood.makeDirty();
        double logLikelihood = likelihood.getLogLikelihood();

        System.out.println("The new and old Likelihood (this value should be roughly the same, debug?): " + logLikelihood + ", " + Double.parseDouble(tree.getAttribute(LIKELIHOOD).toString()));
        if(Double.parseDouble(tree.getAttribute(LIKELIHOOD).toString()) != logLikelihood) {
            /* Newly written check, not sure if this is correct. May need to round values at least */
            //throw new RuntimeException("The values of likelihood are not identical");
        }
//		System.err.printf("New logLikelihood = %4.1f\n", logLikelihood);

        flexTree.setAttribute(LIKELIHOOD, logLikelihood);

        TreeTrait ancestralStates = likelihood.getTreeTrait(TAG);
        for (int i = 0; i < treeModel.getNodeCount(); i++) {

            NodeRef node = treeModel.getNode(i);
            String sample = ancestralStates.getTraitString(treeModel, node);

            String oldSeq = (String) flexTree.getNodeAttribute(flexTree.getNode(i), SEQ_STRING);

            if (oldSeq != null) {

                char[] seq = (sample.substring(1, sample.length() - 1)).toCharArray();

                int length = oldSeq.length();
                //System.out.println("length " + length + "\t" + sample.length() + "\t" + sample + "\t" + seq.length + "\t" + oldSeq);

                for (int j = 0; j < length; j++) {
                    if (oldSeq.charAt(j) == GAP)
                        seq[j] = GAP;
                }

                String newSeq = new String(seq);

//				if( newSeq.contains("MMMMMMM") ) {
//					System.err.println("bad = "+newSeq);
//					System.exit(-1);
//				}

                finalTree.setNodeAttribute(finalTree.getNode(i), NEW_SEQ, newSeq);

            }

//			Taxon taxon = finalTree.getNodeTaxon(finalTree.getNode(i));
//			System.err.println("node: "+(taxon == null ? "null" : taxon.getId())+" "+
//					finalTree.getNodeAttribute(finalTree.getNode(i),NEW_SEQ));

        }
        return finalTree;
    }


    private void setupAttributes(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            Iterator iter = tree.getNodeAttributeNames(node);
            if (iter != null) {
                while (iter.hasNext()) {
                    String name = (String) iter.next();
                    if (!attributeNames.contains(name))
                        attributeNames.add(name);
                }
            }
        }
    }

    private void setupTreeAttributes(Tree tree) {
        Iterator<String> iter = tree.getAttributeNames();
        if (iter != null) {
            while (iter.hasNext()) {
                String name = iter.next();
                treeAttributeNames.add(name);
            }
        }
    }

    private void addTreeAttributes(Tree tree) {
        if (treeAttributeNames != null) {
            if (treeAttributeLists == null) {
                treeAttributeLists = new List[treeAttributeNames.size()];
                for (int i = 0; i < treeAttributeNames.size(); i++) {
                    treeAttributeLists[i] = new ArrayList();
                }
            }

            for (int i = 0; i < treeAttributeNames.size(); i++) {
                String attributeName = treeAttributeNames.get(i);
                Object value = tree.getAttribute(attributeName);

                if (value != null) {
                    treeAttributeLists[i].add(value);
                }
            }
        }

    }


    private Tree summarizeTrees(int burnin, CladeSystem cladeSystem, String inputFileName, boolean useSumCladeCredibility) throws IOException {

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

//		System.out.println("Analyzing " + totalTreesUsed + " trees...");
//		System.out.println("0              25             50             75            100");
//		System.out.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int counter = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName));
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    double score = scoreTree(tree, cladeSystem, useSumCladeCredibility);
//                    System.out.println(score);
                    if (score > bestScore) {
                        bestTree = tree;
                        bestScore = score;
                    }
                }
                if (counter > 0 && counter % stepSize == 0) {
//					System.out.print("*");
//					System.out.flush();
                }
                counter++;
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
//		System.out.println();
//		System.out.println();
        if (useSumCladeCredibility) {
            System.out.println("\tHighest Sum Clade Credibility: " + bestScore);
        } else {
            System.out.println("\tHighest Log Clade Credibility: " + bestScore);
        }

        return bestTree;
    }

    private double scoreTree(Tree tree, CladeSystem cladeSystem, boolean useSumCladeCredibility) {
        if (useSumCladeCredibility) {
            return cladeSystem.getSumCladeCredibility(tree, tree.getRoot(), null);
        } else {
            return cladeSystem.getLogCladeCredibility(tree, tree.getRoot(), null);
        }
    }

    private class CladeSystem {
        //
        // Public stuff
        //

        /**
         */
        public CladeSystem() {
        }

        /**
         * adds all the clades in the tree
         */
        public void add(Tree tree) {
            if (taxonList == null) {
                taxonList = tree;
            }

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is added too (for
            // annotation purposes).
            addClades(tree, tree.getRoot(), null);
            addTreeAttributes(tree);
        }

        public Map<BitSet, Clade> getCladeMap() {
            return cladeMap;
        }

        public Clade getClade(NodeRef node) {
            return null;
        }


        private void addClades(Tree tree, NodeRef node, BitSet bits) {

            BitSet bits2 = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits2.set(index);

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    addClades(tree, node1, bits2);
                }
            }

            addClade(bits2, tree, node);

            if (bits != null) {
                bits.or(bits2);
            }
        }

        private void addClade(BitSet bits, Tree tree, NodeRef node) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(bits);
                cladeMap.put(bits, clade);
            }
            clade.setCount(clade.getCount() + 1);

            if (attributeNames != null) {
                if (clade.attributeLists == null) {
                    clade.attributeLists = new List[attributeNames.size()];
                    for (int i = 0; i < attributeNames.size(); i++) {
                        clade.attributeLists[i] = new ArrayList();
                    }
                }

                for (int i = 0; i < attributeNames.size(); i++) {
                    String attributeName = attributeNames.get(i);
                    Object value;
                    if (attributeName.equals("height")) {
                        value = tree.getNodeHeight(node);
                    } else if (attributeName.equals("length")) {
                        value = tree.getBranchLength(node);
                    } else {
                        value = tree.getNodeAttribute(node, attributeName);
                    }

                    if (value != null) {
                        clade.attributeLists[i].add(value);
                    }
                }
            }
        }

        public void calculateCladeCredibilities(int totalTreesUsed) {
            for (Clade clade : cladeMap.values()) {
                clade.setCredibility(((double) clade.getCount()) / totalTreesUsed);
            }
        }

        public double getSumCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

            double sum = 0.0;

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    sum += getSumCladeCredibility(tree, node1, bits2);
                }

                sum += getCladeCredibility(bits2);

                if (bits != null) {
                    bits.or(bits2);
                }
            }

            return sum;
        }

        public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

            double logCladeCredibility = 0.0;

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
                }

                logCladeCredibility += Math.log(getCladeCredibility(bits2));

                if (bits != null) {
                    bits.or(bits2);
                }
            }

            return logCladeCredibility;
        }

        private double getCladeCredibility(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                return 0.0;
            }
            return clade.getCredibility();
        }

        public void annotateTree(MutableTree tree, NodeRef node, BitSet bits, int heightsOption) {

            BitSet bits2 = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits2.set(index);

                annotateNode(tree, node, bits2, true, heightsOption);
            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    annotateTree(tree, node1, bits2, heightsOption);
                }

                annotateNode(tree, node, bits2, false, heightsOption);
            }

            if (bits != null) {
                bits.or(bits2);
            }
        }

        private void annotateNode(MutableTree tree, NodeRef node, BitSet bits, boolean isTip, int heightsOption) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                throw new RuntimeException("Clade missing");
            }

//			System.err.println("annotateNode new: "+bits.toString()+" size = "+attributeNames.size());
//			for(String string : attributeNames) {
//				System.err.println(string);
//			}
//			System.exit(-1);


            boolean filter = false;
            if (!isTip) {
                double posterior = clade.getCredibility();
                tree.setNodeAttribute(node, "posterior", posterior);
                if (posterior < posteriorLimit) {
                    filter = true;
                }
            }

            for (int i = 0; i < attributeNames.size(); i++) {
                String attributeName = attributeNames.get(i);

                double[] values = new double[clade.attributeLists[i].size()];
                HashMap<String, Integer> hashMap = new HashMap<String, Integer>(clade.attributeLists[i].size());

                if (values.length > 0) {
                    Object v = clade.attributeLists[i].get(0);

                    boolean isHeight = attributeName.equals("height");
                    boolean isBoolean = v instanceof Boolean;

                    boolean isDiscrete = v instanceof String;

                    double minValue = Double.MAX_VALUE;
                    double maxValue = -Double.MAX_VALUE;
                    for (int j = 0; j < clade.attributeLists[i].size(); j++) {
                        if (isDiscrete) {
                            String value = (String) clade.attributeLists[i].get(j);
                            if (value.startsWith("\"")) {
                                value = value.replaceAll("\"", "");
                            }
                            if (attributeName.equals(NEW_SEQ)) {
                                // Strip out gaps before storing
                                value = value.replaceAll("-", "");
                            }
                            if (hashMap.containsKey(value)) {
                                int count = hashMap.get(value);
                                hashMap.put(value, count + 1);
                            } else {
                                hashMap.put(value, 1);

                            }
                        } else if (isBoolean) {
                            values[j] = (((Boolean) clade.attributeLists[i].get(j)) ? 1.0 : 0.0);
                        } else {
                            values[j] = ((Number) clade.attributeLists[i].get(j)).doubleValue();
                            if (values[j] < minValue) minValue = values[j];
                            if (values[j] > maxValue) maxValue = values[j];
                        }
                    }
                    if (isHeight) {
                        if (heightsOption == MEAN_HEIGHTS) {
                            double mean = DiscreteStatistics.mean(values);
                            tree.setNodeHeight(node, mean);
                        } else if (heightsOption == MEDIAN_HEIGHTS) {
                            double median = DiscreteStatistics.median(values);
                            tree.setNodeHeight(node, median);
                        } else {
                            // keep the existing height
                        }
                    }

                    if (!filter) {
                        if (!isDiscrete)
                            annotateMeanAttribute(tree, node, attributeName, values);
                        else
                            annotateModeAttribute(tree, node, attributeName, hashMap);
//						if( tree.getNodeTaxon(node) != null &&
//									tree.getNodeTaxon(node).getId().compareTo("Calanus") == 0) {
//									System.err.println("size = "+hashMap.keySet().size());
//							System.err.println("count = "+hashMap.get(hashMap.keySet().toArray()[0]));
//								}
//							System.err.println();
//						;
                        if (!isBoolean && minValue < maxValue && !isDiscrete) {
                            // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                            // Likewise if it doesn't vary.
                            annotateMedianAttribute(tree, node, attributeName + "_median", values);
                            annotateHPDAttribute(tree, node, attributeName + "_95%_HPD", 0.95, values);
                            annotateRangeAttribute(tree, node, attributeName + "_range", values);
                        }
                    }
                }
            }
        }

        private void annotateMeanAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double mean = DiscreteStatistics.mean(values);
            tree.setNodeAttribute(node, label, mean);
        }

        private void annotateMedianAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double median = DiscreteStatistics.median(values);
            tree.setNodeAttribute(node, label, median);

        }


        public static final String fileName = "junk";
//		public static final String execName = "/sw/bin/clustalw";

        private int externalCalls = 0;

        private void annotateModeAttribute(MutableTree tree, NodeRef node, String label, HashMap<String, Integer> values) {
            if (label.equals(NEW_SEQ)) {
                String consensusSeq = null;
                double[] support = null;
                int numElements = values.keySet().size();
//				System.err.println("size = "+numElements);
                if (numElements > 1) {  // Essentially just finding the modal character
                    try {

                        PrintWriter pw = new PrintWriter(new PrintStream(new FileOutputStream(fileName)));
                        int i = 0;

                        int[] weight = new int[numElements];

//						Iterator iter = values.keySet().iterator();

                        for (String key : values.keySet()) {
                            int thisCount = values.get(key);
                            weight[i] = thisCount;
                            //TODO I THINK I FIXED IT?!?!
                            pw.write(">" + i+"\n");// + " " + thisCount + "\n");
                            pw.write(key + "\n");
                            i++;
                        }
                        pw.close();

                        //Process p = Runtime.getRuntime().exec(kalignExecutable + " " + fileName + " -OUTPUT=NEXUS");
                        //Process p = Runtime.getRuntime().exec(kalignExecutable + " " + fileName + " -fmsf -o"+fileName + ".fasta");
//                        System.out.println("Command: " + kalignExecutable + " " + fileName + " -ffasta -o"+fileName + ".fasta");
//                        Process p = Runtime.getRuntime().exec(kalignExecutable + " " + fileName + " -ffasta -o"+fileName + ".fasta");
                        Process p = Runtime.getRuntime().exec(kalignExecutable + " " + fileName + " -ffasta -q -o"+fileName + ".fasta");

//                        ByteArrayOutputStream kalignOutput = new ByteArrayOutputStream();
//                        StreamGobbler errorGobbler = new StreamGobbler(p.getErrorStream(), "ERR");
//                        StreamGobbler outputGobbler = new StreamGobbler(p.getInputStream(), "OUT", kalignOutput);
//
//                        errorGobbler.start();
//                        outputGobbler.start();
//
//                        int exitVal = p.waitFor();

                        BufferedReader input =
                                new BufferedReader
                                        (new InputStreamReader(p.getInputStream()));
                        {

//                            String line;
                            while ((/*line = */input.readLine()) != null) {
//							System.out.println(line);
                            }
                        }

                        input.close();
                        externalCalls++;
//						System.err.println("clustal call #" + externalCalls);

                        //NexusImporter importer = new NexusImporter(new FileReader(fileName + ".nxs"));
//                        FastaImporter importer = new FastaImporter(new FileReader(new File(fileName + ".fasta")), Nucleotides.INSTANCE);
                        FastaImporter importer = new FastaImporter(new FileReader(new File(fileName + ".fasta")), new GeneralDataType(new String[0]));//AminoAcids.INSTANCE);
                                //new FastaImporter(new FileReader(), new AminoAcids.INSTANCE);

                        Alignment alignment = importer.importAlignment();
                        // build index
                        int[] index = new int[numElements];
                        for (int j = 0; j < numElements; j++)
                            index[j] = alignment.getTaxonIndex("" + j);

                        StringBuffer sb = new StringBuffer();
                        support = new double[alignment.getPatternCount()];

//						System.err.println(new dr.math.matrixAlgebra.Vector(weight));

                        for (int j = 0; j < alignment.getPatternCount(); j++) {
                            int[] pattern = alignment.getPattern(j);
//							support[j] = appendNextConsensusCharacter(alignment,j,sb);
//							System.err.println(new dr.math.matrixAlgebra.Vector(pattern));

                            int[] siteWeight = new int[30]; // + 2 to handle ambiguous and gap
                            int maxWeight = -1;
                            int totalWeight = 0;
                            int maxChar = -1;
                            for (int k = 0; k < pattern.length; k++) {
                                int whichChar = pattern[k];
                                if (whichChar < 30) {
                                    int addWeight = weight[index[k]];
                                    //System.out.println(k + "\t" + addWeight + "\t" + whichChar + "\t" + siteWeight[whichChar]);
                                    siteWeight[whichChar] += addWeight;
                                    totalWeight += addWeight;
//									if( k >= alignment.getDataType().getStateCount()+2 ) {
//										System.err.println("k = "+k);
//										System.err.println("pattern = "+new dr.math.matrixAlgebra.Vector(pattern));
//									}
                                    if (siteWeight[whichChar] > maxWeight) {

                                        maxWeight = siteWeight[whichChar];
                                        maxChar = whichChar;
                                    }
                                } else {
                                    System.err.println("BUG");
                                    System.err.println("k         = " + k);
                                    System.err.println("whichChar = " + whichChar);
                                    System.err.println("pattern   = " + new dr.math.matrixAlgebra.Vector(pattern));
                                }
                            }
                            sb.append(alignment.getDataType().getChar(maxChar));
                            support[j] = (double) maxWeight / (double) totalWeight;
                        }

//						System.err.println("cSeq = " + sb.toString());
                        consensusSeq = sb.toString();

//						System.exit(-1);
                    } catch (FileNotFoundException e) {
                        System.err.println(e.getMessage());
                        System.exit(-1);
                    } catch (Importer.ImportException e) {
                        System.err.println(e.getMessage());
                        System.exit(-1);
                    } catch (IOException e) {
                        System.err.println(e.getMessage());
                        System.exit(-1);
                    }
                } else {
                    consensusSeq = (String) values.keySet().toArray()[0];
                    support = new double[consensusSeq.length()];
                    for (int i = 0; i < support.length; i++)
                        support[i] = 1.0;
                }

                // Trim out gaps from consensus and support
//				ArrayList<Double> newSupport = new ArrayList<Double>(support.length);
                boolean noComma = true;
                StringBuffer newSupport = new StringBuffer("{");
                StringBuffer newSeq = new StringBuffer();
                if (consensusSeq.length() != support.length) {
                    System.err.println("What happened here?");
                    System.exit(-1);
                }

                // Restripping all the sequences of gap characters since they may have been added in during clustal step
                for (int i = 0; i < support.length; i++) {
                    if (consensusSeq.charAt(i) != GAP) {
                        newSeq.append(consensusSeq.charAt(i));
                        if (noComma)
                            noComma = false;
                        else
                            newSupport.append(",");
                        newSupport.append(String.format("%1.3f", support[i]));
                    }
                }
                newSupport.append("}");

                tree.setNodeAttribute(node, label, newSeq);
//				String num = Str
                tree.setNodeAttribute(node, label + ".prob", newSupport);
            } else {
                String mode = null;
                int maxCount = 0;
                int totalCount = 0;

                for (String key : (String[]) values.keySet().toArray()) {
                    int thisCount = values.get(key);
                    if (thisCount == maxCount)
                        mode = mode.concat("+" + key);
                    else if (thisCount > maxCount) {
                        mode = key;
                        maxCount = thisCount;
                    }
                    totalCount += thisCount;
                }
                double freq = (double) maxCount / (double) totalCount;
                tree.setNodeAttribute(node, label, mode);
                tree.setNodeAttribute(node, label + ".prob", freq);
            }
        }

//		private double appendNextConsensusCharacter(Alignment alignment, int j, StringBuffer sb, int[] weight) {
//			int[] pattern = alignment.getPattern(j);
//			int[] siteWeight = new int[alignment.getDataType().getStateCount()];
//			int maxWeight = -1;
//			int totalWeight = 0;
//			int maxChar = -1;
//			for (int k = 0; k < pattern.length; k++) {
//				int whichChar = pattern[k];
//				if (whichChar < alignment.getDataType().getStateCount()) {
//
//					int addWeight = weight[index[k]];
//					siteWeight[whichChar] += addWeight;
//					totalWeight += addWeight;
//					if (siteWeight[k] > maxWeight) {
//						maxWeight = siteWeight[k];
//						maxChar = k;
//					}
//				}
//			}
//			sb.append(alignment.getDataType().getChar(maxChar));
//			return (double) maxWeight / (double) totalWeight;
//
//		}

        private void annotateRangeAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double min = DiscreteStatistics.min(values);
            double max = DiscreteStatistics.max(values);
            tree.setNodeAttribute(node, label, new Object[]{min, max});
        }

        private void annotateHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[] values) {
            int[] indices = new int[values.length];
            HeapSort.sort(values, indices);

            double minRange = Double.MAX_VALUE;
            int hpdIndex = 0;

            int diff = (int) Math.round(hpd * (double) values.length);
            for (int i = 0; i <= (values.length - diff); i++) {
                double minValue = values[indices[i]];
                double maxValue = values[indices[i + diff - 1]];
                double range = Math.abs(maxValue - minValue);
                if (range < minRange) {
                    minRange = range;
                    hpdIndex = i;
                }
            }
            double lower = values[indices[hpdIndex]];
            double upper = values[indices[hpdIndex + diff - 1]];
            tree.setNodeAttribute(node, label, new Object[]{lower, upper});
        }

        class Clade {
            public Clade(BitSet bits) {
                this.bits = bits;
                count = 0;
                credibility = 0.0;
            }

            public int getCount() {
                return count;
            }

            public void setCount(int count) {
                this.count = count;
            }

            public double getCredibility() {
                return credibility;
            }

            public void setCredibility(double credibility) {
                this.credibility = credibility;
            }

            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                final Clade clade = (Clade) o;

                return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

            }

            public int hashCode() {
                return (bits != null ? bits.hashCode() : 0);
            }

            int count;
            double credibility;
            BitSet bits;
            List[] attributeLists = null;
        }

        //
        // Private stuff
        //
        TaxonList taxonList = null;
        Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;
    double posteriorLimit = 0.0;
    String kalignExecutable = null;

    List<String> attributeNames = new ArrayList<String>();
    List<String> treeAttributeNames = new ArrayList<String>();
    List[] treeAttributeLists = null;
    TaxonList taxa = null;

    public static void printTitle() {
        System.out.println();
        centreLine("Ancestral Sequence Annotator " + "v0.1" + ", " + "2008", 60);
//				version.getVersionString() + ", " + version.getDateString(), 60);

        System.out.println();
        centreLine("by", 60);
        System.out.println();
        centreLine("Marc A. Suchard, Wai Lok Sibon Li", 60);
        System.out.println();
        centreLine("Departments of Biomathematics,", 60);
        centreLine("Biostatistics and Human Genetics", 60);
        centreLine("UCLA", 60);
        centreLine("msuchard@ucla.edu", 60);
        System.out.println();
        System.out.println();
        System.out.println("NB: I stole a substantial portion of this code from Andrew Rambaut.");
        System.out.println("    Please also give him due credit.");
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("ancestralsequenceannotator", "<input-file-name> <output-file-name>");
        System.out.println();
        System.out.println("  Example: ancestralsequenceannotator test.trees out.txt");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String targetTreeFileName = null;
        String inputFileName = null;
        String outputFileName = null;


        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        //new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
                        new Arguments.StringOption("heights", new String[]{"keep", "median", "mean"}, false, "an option of 'keep', 'median' or 'mean'"),
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.StringOption("beastInput", new String[]{"true", "false"}, false, "If the input is taken from BEAST rather than BAli-Phy"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
                        new Arguments.Option("help", "option to print this message"),
                        new Arguments.StringOption("kalign", "full_path_to_kalign", "specifies full path to the kalign executable file")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int heights = KEEP_HEIGHTS;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = MEDIAN_HEIGHTS;
            }
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }
        boolean beastInput = false;
        if(arguments.hasOption("beastInput") && arguments.getStringOption("beastInput").equals("true")) {
            SEQ_STRING = "states";
        }

        int target = MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("target")) {
            target = USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }

        // String kalignExecutable = "/usr/local/bin/clustalw";
        String kalignExecutable = "/usr/local/bin/kalign";
        if (arguments.hasOption("kalign")) {
            kalignExecutable = arguments.getStringOption("kalign");
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 2) {
            targetTreeFileName = null;
            inputFileName = args2[0];
            outputFileName = args2[1];
        } else {

            if (inputFileName == null) {
               // No input file name was given so throw up a dialog box...
                inputFileName = Utils.getLoadFileName("AncestralSequenceAnnotator " + version.getVersionString() + " - Select input file file to analyse");
            }
            if (outputFileName == null) {
                outputFileName = Utils.getSaveFileName("AncestralSequenceAnnotator " + version.getVersionString() + " - Select output file");

            }
        }
        if(inputFileName == null || outputFileName == null) {
            System.err.println("Missing input or output file name");
            printUsage(arguments);
            System.exit(1);

        }

        new AncestralSequenceAnnotator(burnin,
                heights,
                posteriorLimit,
                target,
                targetTreeFileName,
                inputFileName,
                outputFileName,
                kalignExecutable);

        System.exit(0);
    }

}

