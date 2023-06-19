package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.trace.TraceException;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * Created by phil on 22/12/2016.
 */
public class GetNSCountsFromTrees {

    private final static Version version = new BeastVersion();
    public static final String BURNIN = "burnin";
    public static final String totalcN = "N";
    public static final String totalcS = "S";
    public static final String totaluN = "b_u_N";
    public static final String totaluS = "b_u_S";
    public static final String historyN = "all_N";
    public static final String historyS = "all_S";
    public static final String SEP = "\t";
    public static final String BRANCHINFO = "branchInfo";
    public static final String[] falseTrue = {"false", "true"};
    public static final String BRANCHSET = "branchSet";
    //    public static final String CLADETAXA = "cladeTaxa";
    public static final String INCLUDECLADES = "includeClades";
    public static final String CLADESTEM = "cladeStem";
    public static final String EXCLUDECLADESTEM = "excludeCladeStem";
    public static final String BACKBONETAXA = "backboneTaxa";
    public static final String ZEROBRANCHES = "zeroBranches";
    public static final String SUMMARY = "summary";
    public static final String SITESUM = "siteSum";
    public static final String CODONSITELIST = "codonSiteList";
    public static final String MRSD = "mrsd";
    public static final String EXCLUDECLADES = "excludeClades";
    public static final String PREFIX = "prefix";
    public static final String HEIGHT_ABOVE = "heightAbove";
    public static final String HEIGHT_BELOW = "heightBelow";
    public static final String TIME_BEFORE = "timeBefore";
    public static final String TIME_AFTER = "timeAfter";
    public static final String DISCRETE_TRAIT_STATE_SET = "discreteTraitStateSet";
    public static final String DISCRETE_TRAIT_ATTRIBUTE = "discreteTraitAttribute";

    public GetNSCountsFromTrees (int burnin, String inputFileName, String outputFileName, boolean branchInfo,
                                 BranchSet branchSet, List<Set> inclusionSets, boolean cladeStem, boolean zeroBranches,
                                 boolean summary, int sites, double mrsd, List<Set> exclusionSets, boolean excludeCladeStems,
                                 double[] siteList, String prefix, double heightAbove, double heightBelow, String discreteTrait, String[] stateSet) throws IOException {

        File inputFile = new File(inputFileName);
        if (inputFile.isFile()){
            System.out.println("Analysing tree file " + inputFileName +" with a burn-in of "+burnin+" trees...");
        } else {
            progressStream.println("cannot find "+inputFileName);
            System.exit(0);
        }

        this.branchInfo = branchInfo;
        // in case you would like to ignore branches without N or S substitutions (can't remember why I thought this could be useful)
        this.zeroBranches = zeroBranches;
        this.summary = summary;
        this.mrsd = mrsd;
        this.cladeStem = cladeStem;
        this.excludeCladeStems = excludeCladeStems;
        this.discreteTrait = discreteTrait;

        this.sites = sites;
        if (siteList!=null){
            this.sites = (siteList.length)*3;
            progressStream.println("Sites set based on site list provided (dN/dS may not be adequate as uN and uS does not account for site list)");
        }

        if (heightAbove > 0 || heightBelow < Double.MAX_VALUE) {
            timeConstraints = true;
        } else {
            timeConstraints = false;
        }

        if (prefix!=null){
            prefixTotalcN = prefix+"."+totalcN;
            prefixTotalcS = prefix+"."+totalcS;
            prefixTotaluN = prefix+"."+totaluN;
            prefixTotaluS = prefix+"."+totaluS;
            prefixHistoryN = prefix+"."+historyN;
            prefixHistoryS = prefix+"."+historyS;
        } else {
            prefixTotalcN = totalcN;
            prefixTotalcS = totalcS;
            prefixTotaluN = totaluN;
            prefixTotaluS = totaluS;
            prefixHistoryN = historyN;
            prefixHistoryS = historyS;
        }

//        resultsStream = System.out;

        if (outputFileName != null) {
            try {
                resultsStream = new PrintStream(new File(outputFileName));
            } catch (IOException e) {
                progressStream.println("Error opening file: " + outputFileName);
                System.exit(-1);
            }
        } else {
            resultsStream = new PrintStream(new File(inputFileName+".NSout.txt"));
        }

        analyze(inputFile, burnin, branchSet, inclusionSets, exclusionSets, siteList, stateSet, heightAbove, heightBelow);

    }

    private void analyze(File inputFile, int burnin, BranchSet branchSet, List<Set> inclusionSets, List<Set> exclusionSets,
                         double[] siteList, String[] stateSet, double heightAbove, double heightBelow){

        if (summary) {
            resultsStream.print("tree"+SEP+"cN"+SEP+"uN"+SEP+"cS"+SEP+"uS"+SEP+"cNrate"+SEP+"cSrate"+SEP+"dN/dS"+"\n");
        }  else {
            resultsStream.print("tree"+SEP+"branch"+SEP+"N/S"+SEP+"site"+SEP+"height/date"+SEP+"fromState"+SEP+"toState");
            if (branchInfo) {
                resultsStream.print(SEP + "branchLength" + SEP + "branchCN/S" + SEP + "branchUN/S" + "\n");
            } else {
                resultsStream.print("\n");
            }
        }

        int totalTrees = 10000;
        int totalStars = 0;

        System.out.println("Reading and analyzing trees (bar assumes 10,000 trees)...");
        System.out.println("0              25             50             75            100");
        System.out.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int count = 0;
        int treeUsed = 1;

        try {
            FileReader fileReader = new FileReader(inputFile);
            TreeImporter importer = new NexusImporter(fileReader);

            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if(count>=burnin) {
                    getNSCounts(tree, treeUsed, branchSet, inclusionSets, exclusionSets, siteList, stateSet, heightAbove, heightBelow);
                    treeUsed ++;
                }
                count++;
                if (totalTrees > 0 && totalTrees % stepSize == 0) {
                    System.out.print("*");
                    totalStars++;
                    if (totalStars % 61 == 0)
                        System.out.print("\n");
                    System.out.flush();
                }
                totalTrees++;
            }
            System.out.print("\n");

        } catch (Importer.ImportException e) {
            progressStream.println("Error Parsing Input Tree: " + e.getMessage());
        } catch (IOException e) {
            progressStream.println("Error Parsing Input Tree: " + e.getMessage());
        }
    }

    private void getNSCounts(Tree tree, int treeUsed, BranchSet branchSet, List<Set> inclusionSets,  List<Set> exclusionSets,
                             double[] siteList, String[] stateSet, double heightAbove, double heightBelow){
        int count = 0;
        double cN = 0;
        double uN = 0;
        double cS = 0;
        double uS = 0;
        double length = 0;

        for (int x = 0; x < tree.getNodeCount(); x++) {
            NodeRef node = tree.getNode(x);
            if (!tree.isRoot(node)){
                count ++;
                if (nodeToConsider(tree, node, branchSet, inclusionSets, exclusionSets, stateSet, heightAbove, heightBelow)){
                    double branchLength = tree.getBranchLength(node);
                    if (timeConstraints){
                        double nodeHeight = tree.getNodeHeight(node);
                        double parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
                        double upperHeight = parentNodeHeight;
                        double lowerHeight = nodeHeight;
                        if (heightBelow < Double.MAX_VALUE) {
                             if (parentNodeHeight > heightBelow) {
                                upperHeight = heightBelow;
                             }
                        }
                        if (heightAbove > 0){
                            if(nodeHeight < heightAbove) {
                                lowerHeight = heightAbove;
                            }
                        }
                        if ((upperHeight - lowerHeight)>0){
                            branchLength = upperHeight - lowerHeight;
                        }
                    }
                    length += branchLength;

                    Object totalNObject = tree.getNodeAttribute(node, prefixTotalcN);
                    Object totalSObject = tree.getNodeAttribute(node, prefixTotalcS);
                    //in case N and S would not be annotated
                    if (totalNObject!=null && totalSObject!=null) {
                        //System.out.println("Hallo");
                        double totalN = (Double) totalNObject;
                        double totalS = (Double) totalSObject;
                        double totaluN = (Double) tree.getNodeAttribute(node, prefixTotaluN);
                        double totaluS = (Double) tree.getNodeAttribute(node, prefixTotaluS);
                        if(siteList==null && !timeConstraints){
                            cN += totalN;
                            cS += totalS;
                        }
                        uN += totaluN;
                        uS += totaluS;
                        if (totalN > 0) {
                            Object[] allNObject = (Object[]) tree.getNodeAttribute(node, prefixHistoryN);
                            for (int a = 0; a < allNObject.length; a++) {
                                Object[] singleNObject = (Object[]) allNObject[a];
                                boolean proceedAgain = false;
                                if(siteList==null) {
                                    if (!timeConstraints){
                                        proceedAgain = true;
                                    } else {
                                        boolean timeCompatible = true;
                                        if (heightBelow < Double.MAX_VALUE){
                                            if((Double) singleNObject[1] > heightBelow){
                                                timeCompatible = false;
                                            }
                                        }
                                        if (heightAbove > 0){
                                            if((Double) singleNObject[1] < heightAbove) {
                                                timeCompatible = false;
                                            }
                                        }
                                        if (timeCompatible){
                                            proceedAgain = true;
                                            cN ++;
                                        }
                                    }
                                } else {
                                    if (inSiteList((Integer)singleNObject[0], siteList)) {
                                        if (!timeConstraints){
                                            proceedAgain = true;
                                            cN ++;
                                        } else {
                                            boolean timeCompatible = true;
                                            if (heightBelow < Double.MAX_VALUE){
                                                if((Double) singleNObject[1] > heightBelow){
                                                    timeCompatible = false;
                                                }
                                            }
                                            if (heightAbove > 0){
                                                if((Double) singleNObject[1] < heightAbove) {
                                                    timeCompatible = false;
                                                }
                                            }
                                            if (timeCompatible){
                                                proceedAgain = true;
                                                cN ++;
                                            }
                                        }
                                    }
                                }
                                if(!summary && proceedAgain){
                                    resultsStream.print(treeUsed + SEP + count + SEP + "N" + SEP);
                                    resultsStream.print(singleNObject[0] + SEP);
                                    if (mrsd > 0) {
                                        resultsStream.print((mrsd - (Double) singleNObject[1]) + SEP);
                                    } else {
                                        resultsStream.print(singleNObject[1] + SEP);
                                    }
                                    resultsStream.print(singleNObject[2] + SEP + singleNObject[3] + SEP);
                                    if (branchInfo) {
                                        resultsStream.print(branchLength + SEP + totalN + SEP + totaluN + "\n");
                                    } else {
                                        resultsStream.print("\n");
                                    }
                                }

                            }
                        }
                        if (totalS > 0) {
                            Object[] allSObject = (Object[]) tree.getNodeAttribute(node, prefixHistoryS);
                            for (int a = 0; a < allSObject.length; a++) {
                                Object[] singleSObject = (Object[]) allSObject[a];
                                boolean proceedAgain = false;
                                if(siteList==null) {
                                    if (!timeConstraints){
                                        proceedAgain = true;
                                    } else {
                                        boolean timeCompatible = true;
                                        if (heightBelow < Double.MAX_VALUE){
                                            if((Double) singleSObject[1] > heightBelow){
                                                timeCompatible = false;
                                            }
                                        }
                                        if (heightAbove > 0){
                                            if((Double) singleSObject[1] < heightAbove) {
                                                timeCompatible = false;
                                            }
                                        }
                                        if (timeCompatible){
                                            proceedAgain = true;
                                            cS ++;
                                        }
                                    }
                                } else {
                                    if (inSiteList((Integer)singleSObject[0], siteList)) {
                                        if (!timeConstraints){
                                            proceedAgain = true;
                                            cS ++;
                                        } else {
                                            boolean timeCompatible = true;
                                            if (heightBelow < Double.MAX_VALUE){
                                                if((Double) singleSObject[1] > heightBelow){
                                                    timeCompatible = false;
                                                }
                                            }
                                            if (heightAbove > 0){
                                                if((Double) singleSObject[1] < heightAbove) {
                                                    timeCompatible = false;
                                                }
                                            }
                                            if (timeCompatible){
                                                proceedAgain = true;
                                                cS ++;
                                            }
                                        }
                                    }
                                }
                                if(!summary && proceedAgain){
                                    resultsStream.print(treeUsed + SEP + count + SEP + "S" + SEP);
                                    resultsStream.print(singleSObject[0] + SEP);
                                    if (mrsd > 0){
                                        resultsStream.print((mrsd - (Double)singleSObject[1]) + SEP);
                                    }  else {
                                        resultsStream.print(singleSObject[1] + SEP);
                                    }
                                    resultsStream.print(singleSObject[2] + SEP + singleSObject[3] + SEP);
                                    if (branchInfo) {
                                        resultsStream.print(branchLength + SEP + totalS + SEP + totaluS + "\n");
                                    } else {
                                        resultsStream.print("\n");
                                    }
                                }
                            }
                        }

                        //in case you do not wish to print out zero branch lengths
                        if ((totalN+totalS)==0){
                            if (zeroBranches) {
                                if (!summary){
                                    resultsStream.print(treeUsed + SEP + count + SEP + "NA" + SEP);
                                    resultsStream.print("NA" + SEP + "NA" + SEP + "NA" + SEP + "NA" + SEP);
                                    if (branchInfo) {
                                        resultsStream.print(branchLength + SEP + "NA" + SEP + "NA" + "\n");
                                    } else {
                                        resultsStream.print("\n");
                                    }
                                }
                            } else {
                                //just to be consistent: if branches with zero N and S are ignored, they should also be ignored in the total branch length sum
                                length -= branchLength;
                            }
                        }
                    }  else {
                        System.err.println("No N or S annotations?");
                        System.exit(-1);
                    }
                }
            }
        }
        if (summary){
            resultsStream.print(treeUsed + SEP + cN + SEP + uN + SEP + cS + SEP + uS
                    + SEP + cN/(length*sites) + SEP + cS/(length*sites) + SEP + (cN/uN)/(cS/uS)+"\n");
        }
    }

    private boolean nodeToConsider(Tree tree, NodeRef node, BranchSet branchSet, List<Set> inclusionSets,  List<Set> exclusionSets,
                                   String[] stateSet, double heightAbove, double heightBelow){
        boolean nodeToConsider = false;
        if (branchSet == BranchSet.ALL) {
            nodeToConsider = true;
        } else if (branchSet == BranchSet.EXT && tree.isExternal(node)) {
            nodeToConsider = true;
        } else if (branchSet == BranchSet.INT && !tree.isExternal(node)) {
            nodeToConsider = true;
        } else if (branchSet == BranchSet.BACKBONE) {
            for (Set inclusionSet: inclusionSets){
                if (onBackbone(tree, node, inclusionSet)){
                    nodeToConsider = true;
                }
            }
        } else if (branchSet == BranchSet.CLADE) {
            for (Set inclusionSet: inclusionSets){
                if (inClade(tree, node, inclusionSet, cladeStem)){
                    nodeToConsider = true;
                }
            }
        } else if (branchSet == BranchSet.SINGLEBRANCH) {
            Set inclusionSet = inclusionSets.get(0); //the issue with multiple inclusion sets is caught in the main method
            if (exclusionSets.size()>0){
                System.err.println("exclusion sets are ignored for a single branch");
            }
            nodeToConsider = isMRCAnode(tree, node, inclusionSet);
        }

        //if the node falls in one of the clades we are excluding, do not proceed.
        if (nodeToConsider && exclusionSets.size()>0 && (branchSet != BranchSet.SINGLEBRANCH)){
            for(Set exclusionSet: exclusionSets){
                if(inClade(tree, node, exclusionSet, excludeCladeStems)){
                    nodeToConsider = false;
                    break;
                }
            }
        }

        // making sure it is a branch that maintains a specified state = node and parent state are a specified state.
        if(nodeToConsider && stateSet!=null) {
            nodeToConsider = false;
            String nodeState = ((String)tree.getNodeAttribute(node, discreteTrait)).replaceAll("\"","");
            String parentNodeState = ((String)tree.getNodeAttribute(tree.getParent(node), discreteTrait)).replaceAll("\"","");;
            if (nodeState.equals(parentNodeState)){
//                System.out.println(nodeState);
                for (String state: stateSet){
                    if (state.equalsIgnoreCase(nodeState)){
                        nodeToConsider = true;
                        break;
                    }
                }
            }
        }

        //partially accounting for specified time constraints
        if (nodeToConsider && timeConstraints) {
            double nodeHeight = tree.getNodeHeight(node);
            double parentNodeHeight = tree.getNodeHeight(tree.getParent(node));
            if(parentNodeHeight < heightAbove){
                nodeToConsider = false;
            }
            if(nodeHeight > heightBelow){
                nodeToConsider = false;
            }
        }

        return nodeToConsider;
    }

    private boolean inSiteList(Integer site, double[] siteList){
        boolean returnBoolean = false;
        for (double siteInList : siteList){
            if (site==siteInList){
                returnBoolean = true;
                break;
            }
        }
        return returnBoolean;
    }

    private static Set getTargetSet(String x) {
        Set targetSet = new HashSet();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(x));
            try {
                String line = reader.readLine().trim();
                while (line != null && !line.equals("")) {
                    targetSet.add(line);
                    line = reader.readLine();
                    if (line != null) line = line.trim();
                }
            }
            catch (IOException io) {
                progressStream.println("Error reading " + x);
            }
        }
        catch (FileNotFoundException a) {
            progressStream.println("Error finding " + x);
        }
        return targetSet;
    }

    private static boolean onBackbone(Tree tree, NodeRef node, Set targetSet) {

        if (tree.isExternal(node)) return false;

        Set leafSet = Tree.Utils.getDescendantLeaves(tree, node);
        int size = leafSet.size();

        leafSet.retainAll(targetSet);

        if (leafSet.size() > 0) {

            // if all leaves below are in target then check just above.
            if (leafSet.size() == size) {

                Set superLeafSet = Tree.Utils.getDescendantLeaves(tree, tree.getParent(node));
                superLeafSet.removeAll(targetSet);

                // the branch is on ancestral path if the super tree has some non-targets in it
                return (superLeafSet.size() > 0);

            } else return true;

        } else return false;
    }

    private static boolean inClade(Tree tree, NodeRef node, Set targetSet, boolean includeStem) {
        Set leafSet = Tree.Utils.getDescendantLeaves(tree, node);

        leafSet.removeAll(targetSet);

        if (leafSet.size() == 0){

            if (includeStem){
                return true;
            }  else {
                Set parentLeafSet = Tree.Utils.getDescendantLeaves(tree, tree.getParent(node));
                parentLeafSet.removeAll(targetSet);
                if (parentLeafSet.size() == 0){
                    return true;
                } else {
                    return false;
                }
            }
        } else {
            return false;
        }
    }

    private static boolean isMRCAnode(Tree tree, NodeRef node, Set targetSet) {
        NodeRef mrca = Tree.Utils.getCommonAncestorNode(tree, targetSet);
        if (node.equals(mrca)){
            return true;
        } else {
            return false;
        }
    }

    private boolean branchInfo;
    private boolean zeroBranches;
    private boolean summary;
    private int sites;
    private double mrsd;
    private boolean cladeStem;
    private boolean excludeCladeStems;
    private String discreteTrait;
    private boolean timeConstraints;
    private static PrintStream progressStream = System.err;
    private PrintStream resultsStream;

    public String prefixTotalcN;
    public String prefixTotalcS;
    public String prefixTotaluN;
    public String prefixTotaluS;
    public String prefixHistoryN;
    public String prefixHistoryS;

    enum BranchSet {
        ALL,
        INT,
        EXT,
        BACKBONE,
        CLADE,
        SINGLEBRANCH
    }

    private static String[] parseVariableLengthStringArray(String inString) {

        List<String> returnList = new ArrayList<String>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            returnList.add(st.nextToken());
        }

        if (returnList.size() > 0) {
            String[] stringArray = new String[returnList.size()];
            stringArray = returnList.toArray(stringArray);
            return stringArray;
        }
        return null;
    }

    public static double[] parseVariableLengthDoubleArray(String inString) {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            returnList.add(Double.parseDouble(st.nextToken()));
        }

        if (returnList.size() > 0) {
            double[] doubleArray = new double[returnList.size()];
            for (int i = 0; i < doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);

            return doubleArray;
        }
        return null;
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
        //arguments.printUsage("loganalyser", "[-burnin <burnin>] [-short][-hpd] [-std] [<input-file-name> [<output-file-name>]]");
        arguments.printUsage("GetNSCountsFromTrees", "[-burnin <burnin>][<input-file-name> [<output-file-name>]]");
        progressStream.println();
        progressStream.println("  Example: GetNSCountsFromTrees -burnin 10000 trees.log out.txt");
        progressStream.println();

    }

    public static void printTitle() {
        System.out.println();
        centreLine("GetNSCountsFromTrees " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey and Marc Suchard", 60);
        System.out.println();
        centreLine("Department of Immunology and Microbiology", 60);
        centreLine("KU Leuven -- University of Leuven", 60);
        centreLine("philippe.lemey@kuleuven.be", 60);
        System.out.println();
        centreLine("Department of Biomathematics", 60);
        centreLine("University of Califormia, Los Angeles", 60);
        centreLine("msuchard@ucla.edu", 60);
        System.out.println();
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption(BRANCHINFO, falseTrue, false, "include a summary for the root [default=off]"),
                        new Arguments.StringOption(BRANCHSET, TimeSlicer.enumNamesToStringArray(BranchSet.values()), false,
                                "branch set [default = all]"),
//                        new Arguments.StringOption(CLADETAXA, "clade taxa file", "specifies a file with taxa that define the clade"),
                        new Arguments.StringOption(INCLUDECLADES, "clade exclusion files", "specifies files with taxa that define clades to be excluded"),
                        new Arguments.StringOption(CLADESTEM, falseTrue, false, "include clade stem [default=false]"),
                        new Arguments.StringOption(EXCLUDECLADESTEM, falseTrue, true, "include clade stem in the exclusion [default=true]"),
                        new Arguments.StringOption(BACKBONETAXA, "Backbone taxa file", "specifies a file with taxa that define the backbone"),
                        new Arguments.StringOption(ZEROBRANCHES, falseTrue, true, "include branches with 0 N and S subtitutions [default=included]"),
                        new Arguments.StringOption(SUMMARY, falseTrue, true, "provide a summary of the N and S counts per tree [default=detailed output]"),
                        new Arguments.RealOption(MRSD, "specifies the most recent sampling data in fractional years to rescale time [default=0]"),
                        new Arguments.StringOption(EXCLUDECLADES, "clade exclusion files", "specifies files with taxa that define clades to be excluded"),
                        new Arguments.IntegerOption(SITESUM, "the number of nucleotide sites to summarize rates in per site per time unit [default = 1]"),
                        new Arguments.StringOption(CODONSITELIST, "list of sites", "sites for which the summary is restricted to"),
                        new Arguments.StringOption(PREFIX, "annotation prefix", "specifies a prefix that is used for the annotations (e.g. to distinguish partition-specific annotations"),

                        new Arguments.RealOption(HEIGHT_ABOVE, "specifies a boundary above which the height must be to be included in the summary [default=none]"),
                        new Arguments.RealOption(HEIGHT_BELOW, "specifies a boundary below which the height must be to be included in the summary [default=none]"),
                        new Arguments.RealOption(TIME_BEFORE, "specifies a boundary before which the time must be to be included in the summary [default=none]"),
                        new Arguments.RealOption(TIME_AFTER, "specifies a boundary after which the time must be to be included in the summary [default=none]"),

                        new Arguments.StringOption(DISCRETE_TRAIT_ATTRIBUTE, "discrete trait attribute", "specifies the string used as attribute for the discrete trait state"),
                        new Arguments.StringOption(DISCRETE_TRAIT_STATE_SET, "discrete trait state set", "specifies which discrete trait states that branches need to maintain to be considered for the summary"),
                        new Arguments.Option("help", "option to print this message")
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

        int burnin = -1;
        if (arguments.hasOption(BURNIN)) {
            burnin = arguments.getIntegerOption(BURNIN);
//            System.err.println("Ignoring a burnin of " + burnin + " trees.");
        }

        boolean branchInfo = false;
        String infoString = arguments.getStringOption(BRANCHINFO);
        if (infoString != null && infoString.compareToIgnoreCase("true") == 0) {
            branchInfo = true;
        }

        BranchSet set = BranchSet.ALL;
        String branch = arguments.getStringOption(BRANCHSET);
        List<Set> inclusionSets = new ArrayList();
        if (branch != null) {
            set = BranchSet.valueOf(branch.toUpperCase());
            progressStream.println("Using the branch set: " + set.name());
        }
        if (set == set.BACKBONE) {
            if (arguments.hasOption(BACKBONETAXA)) {
                String[] fileList = parseVariableLengthStringArray(arguments.getStringOption(BACKBONETAXA));
                for (String singleSet : fileList) {
                    inclusionSets.add(getTargetSet(singleSet));
                    progressStream.println("getting target set for backbone inclusion: " + singleSet);
                }
            } else {
                progressStream.println("you want to get summaries for (a) backbone(s), but no files with taxa to define it are provided??");
            }
        }
        if (set == set.CLADE) {
//            if (arguments.hasOption(CLADETAXA)) {
//                taxaSet = getTargetSet(arguments.getStringOption(CLADETAXA));
//            } else {
//                progressStream.println("you want to get summaries for a clade, but have no taxa to define it??");
//            }

            if (arguments.hasOption(INCLUDECLADES)) {
                String[] fileList = parseVariableLengthStringArray(arguments.getStringOption(INCLUDECLADES));
                for (String singleSet : fileList) {
                    inclusionSets.add(getTargetSet(singleSet));
                    progressStream.println("getting target set for clade inclusion: " + singleSet);
                }
            } else {
                progressStream.println("you want to get summaries for one or more clades, but no files with taxa to define it are provided??");
                System.exit(-1);
            }
        }
        if (set == set.SINGLEBRANCH) {
            if (arguments.hasOption(INCLUDECLADES)) {
                String[] fileList = parseVariableLengthStringArray(arguments.getStringOption(INCLUDECLADES));
                if (fileList.length > 1) {
                    progressStream.println("more than one clade set is specified for a summary of a single branch??");
                    System.exit(-1);
                } else {
                    String singleSet = fileList[0];
                    inclusionSets.add(getTargetSet(singleSet));
                    progressStream.println("getting target set for a single branch summary specification: " + singleSet);
                }
            } else {
                progressStream.println("you want to get summaries for a single branch, but no files with taxa to define to child node of the branch are provided??");
                System.exit(-1);
            }
        }

        boolean cladeStem = false;
        String cladeStemString = arguments.getStringOption(CLADESTEM);
        if (cladeStemString != null && cladeStemString.compareToIgnoreCase("true") == 0) {
            cladeStem = true;
        }

        boolean excludeCladeStems = true;
        String excludeCladeStemString = arguments.getStringOption(CLADESTEM);
        if (excludeCladeStemString != null && excludeCladeStemString.compareToIgnoreCase("false") == 0) {
            excludeCladeStems = false;
        }

        boolean zeroBranches = true;
        String zeroBranchString = arguments.getStringOption(ZEROBRANCHES);
        if (zeroBranchString != null && zeroBranchString.compareToIgnoreCase("false") == 0) {
            zeroBranches = false;
        }

        boolean summary = false;
        String summaryString = arguments.getStringOption(SUMMARY);
        if (summaryString != null && summaryString.compareToIgnoreCase("true") == 0) {
            summary = true;
        }

        int sites = 1;
        if (arguments.hasOption(SITESUM)) {
            sites = arguments.getIntegerOption(SITESUM);
        }

        double mrsd = 0;
        if (arguments.hasOption(MRSD)) {
            mrsd = arguments.getRealOption(MRSD);
        }

        List<Set> exclusionSets = new ArrayList();
        if (arguments.hasOption(EXCLUDECLADES)) {
            String[] fileList = parseVariableLengthStringArray(arguments.getStringOption(EXCLUDECLADES));
            for (String singleSet : fileList) {
                exclusionSets.add(getTargetSet(singleSet));
                progressStream.println("getting target set for clade exclusion: " + singleSet);
            }
        }

        double[] siteList = null;
        if (arguments.hasOption(CODONSITELIST)) {
            siteList = parseVariableLengthDoubleArray(arguments.getStringOption(CODONSITELIST));
            progressStream.println("site list provided: note that dN/dS will not be accurately estimated because the neutral expectation to get dN/dS (uN and uS) is for all sites along a branch.");
        }

        String prefix = arguments.getStringOption(PREFIX);

        double heightAbove = 0;
        if (arguments.hasOption(HEIGHT_ABOVE)) {
            heightAbove = arguments.getRealOption(HEIGHT_ABOVE);
        } else if (arguments.hasOption(TIME_BEFORE)) {
            if (mrsd > 0){
                heightAbove = mrsd - arguments.getRealOption(TIME_BEFORE);
            } else {
                System.err.println(TIME_BEFORE +" is specified but no mrsd (>0)??");
                System.exit(-1);
            }
        }

        double heightBelow = Double.MAX_VALUE;
        if (arguments.hasOption(HEIGHT_BELOW)) {
            heightBelow = arguments.getRealOption(HEIGHT_BELOW);
        } else if (arguments.hasOption(TIME_AFTER)) {
            if (mrsd > 0){
                heightBelow = mrsd - arguments.getRealOption(TIME_AFTER);
            } else {
                System.err.println(TIME_AFTER +" is specified but no mrsd (>0)??");
                System.exit(-1);
            }
        }

        String discreteTrait = null;
        if (arguments.hasOption(DISCRETE_TRAIT_ATTRIBUTE)) {
            discreteTrait = arguments.getStringOption(DISCRETE_TRAIT_ATTRIBUTE);
            progressStream.println("discrete trait attribute provided is "+discreteTrait);
        }

        String[] stateSet = null;
        if (arguments.hasOption(DISCRETE_TRAIT_STATE_SET)) {
            stateSet = parseVariableLengthStringArray(arguments.getStringOption(DISCRETE_TRAIT_STATE_SET));
            if (discreteTrait==null){
                System.err.println("stateSet provided nut no discrete trait attribute provided??");
                System.exit(-1);
            }
            progressStream.println("discrete trait state set provided for summary; first state in set = "+stateSet[0]);
        }

        String inputFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            progressStream.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            inputFileName = args2[0];
        }
        if (args2.length > 1) {
            outputFileName = args2[1];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("GetNSCountsFromTrees " + version.getVersionString() + " - Select tree file to analyse");
        }

        if(burnin==-1) {
            System.err.println("Enter number of trees to burn-in (integer): ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            burnin = Integer.parseInt(br.readLine());
        }

        new GetNSCountsFromTrees(burnin, inputFileName, outputFileName, branchInfo, set, inclusionSets, cladeStem, zeroBranches,
                summary, sites, mrsd, exclusionSets, excludeCladeStems,
                siteList, prefix, heightAbove, heightBelow, discreteTrait, stateSet);

        System.exit(0);
    }

}
