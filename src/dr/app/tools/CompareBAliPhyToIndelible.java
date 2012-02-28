package dr.app.tools;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: sibon
 * Date: 11/16/11
 * Time: 10:16 AM
 * To change this template use File | Settings | File Templates.
 * @author Wai Lok Sibon Li
 */
public class CompareBAliPhyToIndelible {

    public static String CLUSTAL_EXECUTABLE = "/usr/local/bin/clustalw";
    public static final String TEMP_CLUSTAL_FILE_NAME = "clustal_temp";

    public static void main(String args[]) {
        String baliPhyTreeFileName = "";
        String indelibleTreeFileName = "";
        String indelibleAncestralAlignmentFileName = "";

        try {
            baliPhyTreeFileName = args[0];
            indelibleTreeFileName = args[1];
            indelibleAncestralAlignmentFileName = args[2];
            if(args.length > 3) {
                CLUSTAL_EXECUTABLE = args[3];
            }
        } catch(Exception e) {
            System.out.println("Input parameters incorrect");
            System.exit(-1);
        }
        try {
            NexusImporter baliPhyTreeImporter = new NexusImporter(new FileReader(new File(baliPhyTreeFileName)));
            NewickImporter indelibleTreeImporter = new NewickImporter(new FileReader(new File(indelibleTreeFileName)));
            FastaImporter indelibleAlignmentImporter = new FastaImporter(new FileReader(new File(indelibleAncestralAlignmentFileName)), Nucleotides.INSTANCE);


            int baliPhyTreeCount = 0;
            Tree baliPhyTree = new FlexibleTree();
            while (baliPhyTreeImporter.hasTree()) {
                baliPhyTree = baliPhyTreeImporter.importNextTree();
                baliPhyTreeCount++;
            }
            if(baliPhyTreeCount != 1) {
                throw new RuntimeException("There are " + baliPhyTreeCount +
                        " trees in the BAli-Phy tree file. Should only be one. ");
            }

            int indelibleTreeCount = 0;
            Tree indelibleTree = new FlexibleTree();
            while (indelibleTreeImporter.hasTree()) {
                indelibleTree = indelibleTreeImporter.importNextTree();
                indelibleTreeCount++;
            }
            if(indelibleTreeCount != 1) {
                throw new RuntimeException("There are " + indelibleTreeCount +
                        " trees in the BAli-Phy tree file. Should only be one. ");
            }

            if(//indelibleTree.getInternalNodeCount() != baliPhyTree.getInternalNodeCount() ||
                    indelibleTree.getExternalNodeCount() != baliPhyTree.getExternalNodeCount()) {
                throw new RuntimeException("Size of indelible tree not the same as Bali-Phy tree: " +
                //indelibleTree.getInternalNodeCount() + "\t" + baliPhyTree.getInternalNodeCount() + "\t" +
                indelibleTree.getExternalNodeCount() + "\t" + baliPhyTree.getExternalNodeCount() + "\t" +
                indelibleTree.toString() + "\t" + baliPhyTree.toString());
            }

            /* Get the list of taxa under each internal node */
            ArrayList<String>[] taxaListPerNodeBaliPhy = new ArrayList[baliPhyTree.getInternalNodeCount()];
            for(int i=0; i<baliPhyTree.getInternalNodeCount(); i++) {
                Set<String> leaveSet = Tree.Utils.getDescendantLeaves(baliPhyTree, baliPhyTree.getInternalNode(i));
                taxaListPerNodeBaliPhy[i] = new ArrayList(leaveSet);
                Collections.sort(taxaListPerNodeBaliPhy[i]);
            }

            ArrayList<String>[] taxaListPerNodeIndelible = new ArrayList[indelibleTree.getInternalNodeCount()];
            for(int i=0; i<indelibleTree.getInternalNodeCount(); i++) {
                Set<String> leaveSet = Tree.Utils.getDescendantLeaves(indelibleTree, indelibleTree.getInternalNode(i));
                taxaListPerNodeIndelible[i] = new ArrayList(leaveSet);
                Collections.sort(taxaListPerNodeIndelible[i]);
            }


            // TODO
            /* Rooting the unrooted BAli-Phy tree properly */
            //System.out.println("12 days + " + baliPhyTree.getRoot().getNumber());


            /* Getting the matching alignments */
            Alignment indelibleAlignment = indelibleAlignmentImporter.importAlignment();
            int taxaMatchCount = 0;
            for(int i=0; i<taxaListPerNodeBaliPhy.length; i++) {
                for(int j=0; j<taxaListPerNodeIndelible.length; j++) {
                    if(taxaListPerNodeBaliPhy[i].size() == taxaListPerNodeIndelible[j].size()) {
                        boolean taxaMatch = true;
                        for(int k=0; k<taxaListPerNodeIndelible[j].size(); k++) {
                            if(!taxaListPerNodeBaliPhy[i].get(k).equals(taxaListPerNodeIndelible[j].get(k))) {
                                taxaMatch = false;
                            }
                        }

                        if(taxaMatch) {
                            String baliPhyAncestralSequence = (String) baliPhyTree.getNodeAttribute(
                                    baliPhyTree.getInternalNode(i), "newSeq");

                            String indelibleSequenceName = (String) indelibleTree.getNodeAttribute(indelibleTree.getInternalNode(j), "label");
                            System.out.print(indelibleSequenceName + "\t");


                            //System.out.println("d " + indelibleAlignment.getTaxon(1).getId());
                            int taxonIndex = indelibleAlignment.getTaxonIndex(indelibleSequenceName);
                            String gappedIndelibleSequence = indelibleAlignment.getAlignedSequenceString(taxonIndex);

//                            for(int k=0; k<taxaListPerNodeIndelible[j].size(); k++) {
//                                System.out.println("|" + taxaListPerNodeBaliPhy[i].get(k) + "\t" + (taxaListPerNodeIndelible[j].get(k)));
//                            }

                            String ungappedIndelibleSequence = gappedIndelibleSequence.replaceAll("\\*", "").replaceAll("-", "");

                            calculateStatistics(baliPhyAncestralSequence, ungappedIndelibleSequence);

                            //System.out.println(baliPhyAncestralSequence + "\t" + ungappedIndelibleSequence);

                            taxaMatchCount++;
                        }
                    }
                }

            }
            System.out.println("Taxa match count: " + taxaMatchCount + "/" + taxaListPerNodeIndelible.length);


//            String[] taxonSetStrings = sb.toString().split("\t");
//            for(int j=0; j<taxonSetStrings.length; j++) {
//                int taxonIndex = indelibleAlignment.getTaxonIndex(taxonSetStrings[j]);
//                String indelibleSequence = indelibleAlignment.getAlignedSequenceString(taxonIndex);
//
//                String baliPhySequence = taxonSetStrings[j]
//            }


            //BufferedReader br = new BufferedReader(new FileReader(new File(baliPhyTreeFileName)));
            //BufferedReader br = new BufferedReader(new FileReader(new File(indelibleTreeFileName)));
            //BufferedReader br = new BufferedReader(new FileReader(new File(indelibleAncestralAlignmentFileName)));
        }
        catch(IOException e) {
            e.printStackTrace();
        }
        catch(Importer.ImportException e) {
            e.printStackTrace();
        }
    }




    public static void calculateStatistics(String sequence1, String sequence2) {

        // Compute alignment using CLUSTAL, see ASA
        String alignmentSeq1 = "";
        String alignmentSeq2 = "";

        try {

            File f = new File(TEMP_CLUSTAL_FILE_NAME + ".fasta");

            PrintWriter pw = new PrintWriter(new PrintStream(f));


            pw.write(">" + "sequence1" + "\n");
            pw.write(sequence1+"\n");
            pw.write(">" + "sequence2" + "\n");
            pw.write(sequence2+"\n");

            pw.close();



            Process p = Runtime.getRuntime().exec(CLUSTAL_EXECUTABLE + " " + TEMP_CLUSTAL_FILE_NAME + ".fasta" + " -OUTPUT=NEXUS");
            BufferedReader input =
                                new BufferedReader
                                        (new InputStreamReader(p.getInputStream()));
                        {
                            //String line;

                            while ((/*line = */input.readLine()) != null) {
							    //System.out.println(line);
                            }
                        }
            input.close();
            BufferedReader err =
                                new BufferedReader
                                        (new InputStreamReader(p.getErrorStream()));
                        {
                            String line;

                            while ((line = err.readLine()) != null) {
							    System.err.println(line);
                            }
                        }
            err.close();

//            BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
////            {
////                //String line;
////                while ((/*line = */input.readLine()) != null) {
////    //							System.out.println(line);
////                }
////            }
//            input.close();

            NexusImporter importer = new NexusImporter(new FileReader(TEMP_CLUSTAL_FILE_NAME + ".nxs"));

            f.delete();

            Alignment alignment = importer.importAlignment();
            alignmentSeq1 = alignment.getAlignedSequenceString(alignment.getTaxonIndex("sequence1"));
            alignmentSeq2 = alignment.getAlignedSequenceString(alignment.getTaxonIndex("sequence2"));

        } catch (Importer.ImportException e) {
            System.err.println(e.getMessage());
            System.exit(-1);
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(-1);
        }

        char[] alignment1 = alignmentSeq1.toCharArray();
        char[] alignment2 = alignmentSeq2.toCharArray();
        if(alignment1.length != alignment2.length) {
            throw new RuntimeException("Length of sequences from alignment are not the same");
        }

        double sequenceIdentity = calculateSequenceIdentity(alignment1, alignment2);
        System.out.print(sequenceIdentity + "\t");
        int gapsMissed = calculateGapsMissed(alignment1, alignment2);
        System.out.print(gapsMissed + "\t");
        int sequenceLength = calculateSequenceLength(sequence1, sequence2);
        System.out.print(sequenceLength + "\t");
        System.out.println(alignmentSeq1 + "\t" + alignmentSeq2);

    }

    public static int calculateSequenceLength(String sequence1, String sequence2) {
        return sequence1.length() - sequence2.length();

    }

    public static double calculateSequenceIdentity(char[] alignment1, char[] alignment2) {
        // number of characters matching divided by total length of alignment
        int count = 0;
        for(int i=0; i<alignment1.length; i++) {
            if(alignment1[i] == alignment2[i]) {
                //System.out.println(alignment1[i] + "\t" + alignment2[i]);
                count++;
            }
        }
        //System.out.println("hello " + count + "\t" + alignment1.length);
        return count*1.0/alignment1.length;

    }

    public static int calculateGapsMissed(char[] alignment1, char[] alignment2) {
        // Number of gap characters that the character misplaced.
        int count = 0;
        for(int i=0; i<alignment1.length; i++) {
            if((alignment1[i]=='-' && alignment2[i]!='-') || (alignment2[i]=='-' && alignment1[i]!='-')) {
                count++;
            }
        }
        return count;
    }






//    private void getCladeSubset(Tree tree, NodeRef node, ArrayList<String> taxaSet) {
//        if(tree.isExternal(node))  {
//
//            //taxaSet.append( + "\t");
//            String taxaName = tree.getNodeTaxon(node).getId();
//            taxaSet.add(taxaName);
//            //taxaSet.add(tree.getNodeTaxon(node));
//            //taxaSet.append(tree.getTaxonIndex(tree.getNodeTaxon(node).getId()) + "\t");
//            //taxaSet.add(tree.getTaxonIndex(tree.getNodeTaxon(node)).getId());
//        }
//        else {
//            for (int i = 0; i < tree.getChildCount(node); i++) {
//                NodeRef n = tree.getChild(node, i);
//                getCladeSubset(tree, n, taxaSet);
//            }
//        }
//
//    }

}
