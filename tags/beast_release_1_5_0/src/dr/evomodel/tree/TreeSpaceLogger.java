/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author shhn001
 */
public class TreeSpaceLogger {

    private String[] trees;
    private int[] islands;
    private int[][] transitions;
    private int[][] steps;
    private int[][] islandTransitions;
    private int[] islandDwellTimes;
    private String lastKnownTree;
    private List<String> currentPath;
//	private HashMap<String, Integer>[][] pathes;

    /**
     *
     */
    public TreeSpaceLogger(String filename) {
        readTrees(filename);
        currentPath = new ArrayList<String>();
    }

    private void readTrees(String filename) {
        try {
            FileReader fr = new FileReader(new File(filename));
            BufferedReader br = new BufferedReader(fr);
            String line;
            List<Tree> trees = new ArrayList<Tree>();
            List<Integer> islands = new ArrayList<Integer>();

            while ((line = br.readLine()) != null) {
                // read trees
                String[] tokens = line.split("\\s++");

                NewickImporter importer = new NewickImporter(line);
                Tree t = importer.importNextTree();
                trees.add(t);
                islands.add(Integer.parseInt(tokens[0]));
            }

            br.close();
            fr.close();

            this.trees = new String[trees.size()];
            this.islands = new int[trees.size()];
            int maxIsland = 0;
            for (int i = 0; i < trees.size(); i++) {
                Tree tree = trees.get(i);
                String newick = Tree.Utils.uniqueNewick(tree, tree.getRoot());
                this.trees[i] = newick;
                this.islands[i] = islands.get(i) - 1;
                if (islands.get(i) > maxIsland) {
                    maxIsland = islands.get(i);
                }
            }

            transitions = new int[trees.size()][trees.size()];
            steps = new int[trees.size()][trees.size()];
            islandTransitions = new int[maxIsland][maxIsland];
            islandDwellTimes = new int[maxIsland];

//			pathes = new HashMap[trees.size()][trees.size()];
//			for (int i = 0; i < trees.size(); i++) {
//				for (int j = 0; j < trees.size(); j++) {
//					pathes[i][j] = new HashMap<String, Integer>();
//				}
//			}

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (ImportException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void log(Tree tree) {
        String currentTree = Tree.Utils.uniqueNewick(tree, tree.getRoot());

        if (lastKnownTree == null || lastKnownTree.equals("")) {
            // check if we ended up at another known tree
            int index = findTree(currentTree);
            if (index != -1) {
                lastKnownTree = currentTree;
            }
        } else {
            // check if we left our current tree or if went in a unlikely circle
            if (currentTree.equals(lastKnownTree)) {
                if (currentPath.size() >= 1) {
                    int startIndex = findTree(lastKnownTree);
                    transitions[startIndex][startIndex]++;
                    steps[startIndex][startIndex] += currentPath.size();
                    islandTransitions[islands[startIndex]][islands[startIndex]]++;
                    islandDwellTimes[islands[startIndex]] += currentPath.size();
                }
                currentPath.clear();
            } else {
                // check if we ended up at another known tree
                int index = findTree(currentTree);
                if (index == -1) {
                    // we didn't end at a know tree yet
                    currentPath.add(currentTree);
                } else {
                    currentPath.add(currentTree);
                    int startIndex = findTree(lastKnownTree);
                    transitions[startIndex][index]++;
                    steps[startIndex][index] += currentPath.size() - 1;
                    islandDwellTimes[islands[index]] += currentPath.size() - 1;
                    islandTransitions[islands[startIndex]][islands[index]]++;

//					if (pathes[startIndex][index].containsKey(currentPath
//							.toString())) {
//						int trans = pathes[startIndex][index].get(currentPath
//								.toString());
//						pathes[startIndex][index].put(currentPath.toString(),
//								trans + 1);
//					} else {
//						pathes[startIndex][index]
//								.put(currentPath.toString(), 1);
//					}

                    currentPath.clear();
                    lastKnownTree = currentTree;
                }
            }
        }
    }

    private int findTree(String tree) {
        for (int i = 0; i < trees.length; i++) {
            if (trees[i].equals(tree)) {
                return i;
            }
        }
        return -1;
    }

    public void printTransitionMatrix() {
        System.out.println("\n\nPrint transition matrix:\n");

        for (int i = 0; i < trees.length; i++) {
            for (int j = 0; j < trees.length; j++) {
                System.out.print(transitions[i][j] + "\t");
            }
            System.out.println();
        }
    }

//	public void printPathesMatrix() {
//		System.out.println("\n\nPrint # of pathes matrix:\n");
//
//		for (int i = 0; i < trees.length; i++) {
//			for (int j = 0; j < trees.length; j++) {
//				System.out.print(pathes[i][j].size() + "\t");
//			}
//			System.out.println();
//		}
//	}

    public void printIslandTransitionMatrix() {
        System.out.println("\n\nPrint island transition matrix:\n");

        for (int i = 0; i < islandTransitions.length; i++) {
            for (int j = 0; j < islandTransitions.length; j++) {
                System.out.print(islandTransitions[i][j] + "\t");
            }
            System.out.println();
        }
    }

    public void printMeanTopologyMatrix() {
        System.out.println("\n\nPrint mean # of topologies on path matrix:\n");

        for (int i = 0; i < trees.length; i++) {
            for (int j = 0; j < trees.length; j++) {
                int trans = transitions[i][j];
                if (trans == 0) {
                    trans = 1;
                }
                System.out.print((steps[i][j] / trans) + "\t");
            }
            System.out.println();
        }
    }

    public void printIslandDwellMatrix() {
        System.out.println("\n\nPrint mean dwell time matrix:\n");

        for (int i = 0; i < islandDwellTimes.length; i++) {
            int transitions = 0;
            for (int j = 0; j < islandDwellTimes.length; j++) {
                if (i != j) {
                    transitions += islandTransitions[i][j];
                }
            }
            if (transitions == 0) {
                transitions = 1;
            }
            System.out.print((islandDwellTimes[i] / transitions) + "\t");
            System.out.println();
        }
    }

    public void finishLogging() {
        printTransitionMatrix();
//		printPathesMatrix();
        printIslandTransitionMatrix();
        printMeanTopologyMatrix();
        printIslandDwellMatrix();
    }

    public static void main(String[] args) {
        TreeSpaceLogger logger = new TreeSpaceLogger("DS1_Treespace.trees");
        logger.finishLogging();
    }

}
