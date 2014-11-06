package dr.app.tools;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Takes a tree with a complete jump history annotated to its branches and makes new internal nodes of degree 1
 * corresponding to the jump points.
 *
 * @author Matthew Hall
 */
public class BranchJumpPlotter {

    private NexusImporter treesIn;
    private NexusExporter treesOut;
    private String traitName;

    public BranchJumpPlotter(NexusImporter treesIn, NexusExporter treesOut, String traitName){
        this.treesIn = treesIn;
        this.treesOut = treesOut;
        this.traitName = traitName;
    }

    private FlexibleTree rewireTree(Tree tree, boolean verbose){
        int totalJumps = 0;
        FlexibleTree outTree = new FlexibleTree(tree, true);
        for(int nodeNo = 0; nodeNo < outTree.getNodeCount(); nodeNo++){
            FlexibleNode node = (FlexibleNode)outTree.getNode(nodeNo);
            String finalHost = (String)node.getAttribute(traitName);
            node.setAttribute(traitName,finalHost.replaceAll("\"",""));
            HashSet<String[]> jumps = readCJH(node);
            if(verbose){
                System.out.print("Node "+nodeNo+": ");
            }
            if(jumps!=null){
                double[] heights = new double[jumps.size()];
                HashMap<Double,String[]> heightMap = new HashMap<Double, String[]>();
                int count = 0;
                for(String[] jump: jumps){
                    double height = Double.parseDouble(jump[1]);
                    heights[count] = height;
                    heightMap.put(height,jump);
                    count++;
                }
                Arrays.sort(heights);
                FlexibleNode oldParent = (FlexibleNode)tree.getParent(node);
                FlexibleNode needsNewParent = node;

                for (double height : heights) {
                    totalJumps++;
                    if(!needsNewParent.getAttribute(traitName).equals(heightMap.get(height)[3])){
                        throw new RuntimeException("Destination traits do not match");
                    }
                    FlexibleNode parent = (FlexibleNode)outTree.getParent(needsNewParent);
                    outTree.beginTreeEdit();
                    outTree.removeChild(parent, needsNewParent);
                    needsNewParent.setLength(height-needsNewParent.getHeight());
                    FlexibleNode jumpNode = new FlexibleNode();
                    jumpNode.setHeight(height);
                    jumpNode.setLength(parent.getHeight() - height);
                    jumpNode.setAttribute(traitName,heightMap.get(height)[2]);
                    outTree.addChild(parent, jumpNode);
                    outTree.addChild(jumpNode, needsNewParent);
                    outTree.endTreeEdit();
                    needsNewParent = jumpNode;
                }
            }
            if(verbose){
                if(jumps==null){
                    System.out.println(0+" ("+totalJumps+")");
                } else {
                    System.out.println(jumps.size()+" ("+totalJumps+")");
                }
            }
        }

        outTree = new FlexibleTree((FlexibleNode)outTree.getRoot());
        if(verbose){
            System.out.println("Total jumps: "+totalJumps);
            int[] childCounts = new int[3];
            for(int i=0; i<outTree.getNodeCount(); i++){
                childCounts[outTree.getChildCount(outTree.getNode(i))]++;
            }
            for(int i=0; i<3; i++){
                System.out.println(childCounts[i]+" nodes have "+i+" children");
            }
        }
        return outTree;
    }

    private HashSet<String[]> readCJH(FlexibleNode node){
        if(node.getAttribute("history_all")!=null){
            HashSet<String[]> out = new HashSet<String[]>();
            String cjh = (String)node.getAttribute("history_all");
            cjh = cjh.replaceFirst("<<","");
            cjh = cjh.replaceFirst(">>","");
            String[] blocks = cjh.split(">\\|<");
            for(String block: blocks){
                String[] items = block.split("\\|");
                out.add(items);
            }
            return out;
        } else {
            return null;
        }
    }

    private void translateTreeFile(){
        try{
            ArrayList<Tree> trees = new ArrayList<Tree>();
            int count = 1;
            while(treesIn.hasTree()){
                if(count % 100 == 0){
                    System.out.println("Doing tree "+count);
                }
                trees.add(rewireTree(treesIn.importNextTree(),false));
                count++;
            }
            Tree[] treeArray = trees.toArray(new Tree[trees.size()]);
            treesOut.exportTrees(treeArray);
        } catch(IOException e){
            System.out.println("Problem reading file ("+e.toString()+")");
        } catch(Importer.ImportException e){
            System.out.println("Problem importing trees ("+e.toString()+")");
        }
    }

    public static void main(String[] args){
        try{
            String traitName = args[0];
            String inputFile = args[1];
            String outputFile = args[2];
            NexusImporter importer = new NexusImporter(new FileReader(inputFile));
            NexusExporter exporter = new NexusExporter(new PrintStream(outputFile));
            BranchJumpPlotter plotter = new BranchJumpPlotter(importer,exporter,traitName);
            plotter.translateTreeFile();
        } catch(FileNotFoundException e){
            System.out.println("File not found");
        }

    }



}
