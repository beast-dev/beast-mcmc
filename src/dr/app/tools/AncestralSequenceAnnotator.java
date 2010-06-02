package dr.app.tools;

import dr.app.util.Arguments;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeExporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.WAG;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.inference.model.Parameter;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/*
 * @author Marc A. Suchard
 */

public class AncestralSequenceAnnotator {

    public final static int MAX_CLADE_CREDIBILITY = 0;
    public final static int MAX_SUM_CLADE_CREDIBILITY = 1;
    public final static int USER_TARGET_TREE = 2;

    public final static int KEEP_HEIGHTS = 0;
    public final static int MEAN_HEIGHTS = 1;
    public final static int MEDIAN_HEIGHTS = 2;

    public AncestralSequenceAnnotator(int burnin,
                                      int heightsOption,
                                      double posteriorLimit,
                                      int targetOption,
                                      String targetTreeFileName,
                                      String inputFileName,
                                      String outputFileName,
                                      String clustalExecutable
    ) throws IOException {

        this.posteriorLimit = posteriorLimit;

        this.clustalExecutable = clustalExecutable;

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
                    tree = unprocessedTree;
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

//		System.exit(-1);
        System.out.println("Writing annotated tree....");
        if (outputFileName != null) {
            exporter = new NexusExporter(new PrintStream(new FileOutputStream(outputFileName)));
            exporter.exportTree(targetTree);
        } else {
            exporter = new NexusExporter(System.out);
            exporter.exportTree(targetTree);
        }


    }

    abstract class SubstitutionModelLoader {

        protected SubstitutionModel substModel;
        protected FrequencyModel freqModel;
        protected SimpleAlignment alignment;
        private final String name;

        public SubstitutionModel getSubstitutionModel() {
            return substModel;
        }

        public FrequencyModel getFrequencyModel() {
            return freqModel;
        }

        public SimpleAlignment getAlignment() {
            return alignment;
        }

        public String getName() {
            return name;
        }

        protected abstract void modelSpecifics(Tree tree);

        public void load(Tree tree) {
            alignment = new SimpleAlignment();
            modelSpecifics(tree);
            printLogLikelihood(tree);
        }

        SubstitutionModelLoader(String name) {
            this.name = name;
        }

        protected boolean doPrint = false;

        private void printLogLikelihood(Tree tree) {
            if (doPrint) {
                Double logLikelihood = (Double) tree.getAttribute(LIKELIHOOD);
                if (logLikelihood != null)
                    System.err.printf("%5.1f", logLikelihood);
            }
        }

        protected void doPrint(Tree tree) {

        }

    }

    public SubstitutionModelLoader[] modelLoaders = {
            new SubstitutionModelLoader("Empirical(wag.dat)*pi") {
                protected void modelSpecifics(Tree tree) {
                    // Instantiate a WAG model
                    double[] freq = new double[20];
                    int cnt = 0;
                    double sum = 0;
                    String charList = "";
                    for (Iterator<String> i = tree.getAttributeNames(); i.hasNext();) {
                        String name = i.next();
                        if (name.startsWith("pi")) {
                            String character = name.substring(2, 3);
                            charList = charList.concat(character);
                            Double value = (Double) tree.getAttribute(name);
                            freq[cnt++] = value;
                            sum += value;
                        }
                    }
                    double[] freqOrdered = new double[20];
                    for (int i = 0; i < 20; i++) {
                        int index = charList.indexOf(AA_ORDER.charAt(i));
                        freqOrdered[i] = freq[index] / sum;
                    }
                    freqModel = new FrequencyModel(AminoAcids.INSTANCE, new Parameter.Default(freqOrdered));
                    substModel = new EmpiricalAminoAcidModel(WAG.INSTANCE, freqModel);
                    alignment.setDataType(AminoAcids.INSTANCE);
//					doPrint = true;
                }
            }
    };


    public static final String KAPPA_STRING = "kappa";
    public static final String SEQ_STRING = "seq";
    public static final String NEW_SEQ = "newSeq";
    public static final String TAG = "tag";
    public static final String LIKELIHOOD = "lnL";
    public static final String SUBST_MODEL = "subst";
    public static final String AA_ORDER = "ACDEFGHIKLMNPQRSTVWY";

    //	public static final String WAG_STRING = "Empirical(Data/wag.dat)*pi";
    private final int die = 0;

    private Tree processTree(Tree tree) {

        // Remake tree to fix node ordering

        String modelType = (String) tree.getAttribute(SUBST_MODEL);
        SubstitutionModelLoader loader = null;

        for (int i = 0; i < modelLoaders.length && loader == null; i++) {
            if (modelType.equals(modelLoaders[i].getName())) {
                loader = modelLoaders[i];
            }
        }

        if (loader == null) {
            System.err.println("Substitution model type '" + modelType + "' not implemented");
            System.exit(-1);
        }

        loader.load(tree);

        SubstitutionModel substModel = loader.getSubstitutionModel();
        SimpleAlignment alignment = loader.getAlignment();

        // Get sequences
        String[] sequence = new String[tree.getNodeCount()];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            sequence[i] = (String) tree.getNodeAttribute(node, SEQ_STRING);
            if (tree.isExternal(node)) {
                Taxon taxon = tree.getNodeTaxon(node);
                alignment.addSequence(new Sequence(taxon, sequence[i]));
            }
        }

        // Make evolutionary model
        SiteModel siteModel = new GammaSiteModel(substModel, new Parameter.Default(1.0), null, 0, null);
        BranchRateModel rateModel = new StrictClockBranchRates(new Parameter.Default(1.0));

        FlexibleTree flexTree = sampleTree(tree, alignment, siteModel, rateModel);
        introduceGaps(flexTree, tree);

        return flexTree;
    }

    private void introduceGaps(FlexibleTree flexTree, Tree gapTree) {
        // I forget what this function was supposed to do.
    }


    public static final char GAP = '-';

    boolean[] bit = null;

    private FlexibleTree sampleTree(Tree tree, PatternList alignment, SiteModel siteModel, BranchRateModel rateModel) {
        FlexibleTree flexTree = new FlexibleTree(tree, true);
        flexTree.adoptTreeModelOrdering();
        FlexibleTree finalTree = new FlexibleTree(tree);
        finalTree.adoptTreeModelOrdering();
        TreeModel treeModel = new TreeModel(tree);

        // Turn off noisy logging by TreeLikelihood constructor
        Logger logger = Logger.getLogger("dr.evomodel");
        boolean useParentHandlers = logger.getUseParentHandlers();
        logger.setUseParentHandlers(false);

        AncestralStateTreeLikelihood likelihood = new AncestralStateTreeLikelihood(
                alignment,
                treeModel,
                siteModel,
                rateModel,
                false, true,
                alignment.getDataType(),
                TAG,
                false);

        logger.setUseParentHandlers(useParentHandlers);

        // Sample internal nodes
        likelihood.makeDirty();
        double logLikelihood = likelihood.getLogLikelihood();

//		System.err.printf("New logLikelihood = %4.1f\n", logLikelihood);

        flexTree.setAttribute(LIKELIHOOD, logLikelihood);

        TreeTrait ancestralStates = likelihood.getTreeTrait(AncestralStateTreeLikelihood.STATES_KEY);
        for (int i = 0; i < treeModel.getNodeCount(); i++) {

            NodeRef node = treeModel.getNode(i);
            String[] sample = ancestralStates.getTraitString(treeModel, node);

            String oldSeq = (String) flexTree.getNodeAttribute(flexTree.getNode(i), SEQ_STRING);

            if (oldSeq != null) {

                char[] seq = (sample[0].substring(1, sample[0].length() - 1)).toCharArray();

                int length = oldSeq.length();

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
//					System.err.println("attr = "+attributeName);
                    Object value;
                    if (attributeName.equals("height")) {
                        value = tree.getNodeHeight(node);
                    } else if (attributeName.equals("length")) {
                        value = tree.getBranchLength(node);
                    } else {
                        value = tree.getNodeAttribute(node, attributeName);
//						System.err.println(attributeName+" : "+(String)value);
                    }

                    //if (value == null) {
                    //    System.out.println("attribute " + attributeNames[i] + " is null.");
                    //}

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
//				System.err.println(attributeName);

                double[] values = new double[clade.attributeLists[i].size()];
                HashMap<String, Integer> hashMap = new HashMap<String, Integer>(clade.attributeLists[i].size());
//				Hashtable hashMap = new Hashtable(clade.attributeLists[i].size());


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
//								System.err.println("good");
                            } else {
                                hashMap.put(value, 1);

                            }
//							if( tree.getNodeTaxon(node) != null &&
//									tree.getNodeTaxon(node).getId().compareTo("Calanus") == 0) {
//									System.err.println(value);
//								}
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
//				if( tree.isExternal(node) ) {
//
//					return;
//				}
                String consensusSeq = null;
                double[] support = null;
                int numElements = values.keySet().size();
//				System.err.println("size = "+numElements);
                if (numElements > 1) {
//					System.err.println("More than one!");
                    try {

                        PrintWriter pw = new PrintWriter(new PrintStream(new FileOutputStream(fileName)));
                        int i = 0;

                        int[] weight = new int[numElements];

//						Iterator iter = values.keySet().iterator();

                        for (String key : values.keySet()) {
                            int thisCount = values.get(key);
                            weight[i] = thisCount;
                            pw.write(">" + i + " " + thisCount + "\n");
                            pw.write(key + "\n");
                            i++;
                        }
                        pw.close();

                        Process p = Runtime.getRuntime().exec(clustalExecutable + " " + fileName + " -OUTPUT=NEXUS");
                        BufferedReader input =
                                new BufferedReader
                                        (new InputStreamReader(p.getInputStream()));
                        {
                            //String line;

                            while ((/*line = */input.readLine()) != null) {
//							System.out.println(line);
                            }
                        }

                        input.close();
                        externalCalls++;
//						System.err.println("clustal call #" + externalCalls);

                        NexusImporter importer = new NexusImporter(new FileReader(fileName + ".nxs"));

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
    String clustalExecutable = null;

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
        centreLine("Marc A. Suchard", 60);
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
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annoated"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
                        new Arguments.Option("help", "option to print this message"),
                        new Arguments.StringOption("clustal", "full_path_to_clustal", "specifies full path to the clustal executable file")
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

        int target = MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("target")) {
            target = USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }

        String clustalExecutable = "/usr/local/bin/clustalw";
        if (arguments.hasOption("clustal")) {
            clustalExecutable = arguments.getStringOption("clustal");
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
                clustalExecutable);

        System.exit(0);
    }

}

