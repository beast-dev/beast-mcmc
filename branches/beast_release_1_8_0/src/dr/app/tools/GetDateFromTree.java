package dr.app.tools;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.Tree;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * @author Walter Xie
 */
public class GetDateFromTree extends NewickImporter {

    public GetDateFromTree(Reader reader) {
        super(reader);
    }

    public GetDateFromTree(String treeString) {
        super(treeString);
    }

    static public void main(String[] args) {
        for (int c = 1; c <= 10; c++) {
            int index = 10;

            try {
                System.out.println("Input trees from " + pathInput + c + inputFileName + "\n\n");
                FileReader fileReader = new FileReader(pathInput + c + inputFileName);
                GetDateFromTree getDateFromTree = new GetDateFromTree(fileReader); // many trees
                int totalTrees = 0;

                try {
                    while (getDateFromTree.hasTree()) {
                        Tree tree = getDateFromTree.importNextTree();
//                        System.out.println(tree.toString());
//                        System.out.println(insetTreeIndex(2, tree.toString()));

                        if (totalTrees == index * 100) { // 1000
                            System.out.println("input " + totalTrees + "th tree");
                            getDate(Integer.toString(c), index, tree);
                            index += 10;
                        }
                        totalTrees++;
                    }

                } catch (ImportException e) {
                    System.err.println("Error Parsing Input Tree: " + e.getMessage());
                    return;
                }

                fileReader.close();
                System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++\n\n");

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
//        for (String t : trees) {
//            GetDateFromTree getDateFromTree = new GetDateFromTree(t);
//            getDate(index, getDateFromTree, t);
//            index += 10;
//        }
    }

    private static void getDate(String curD, int index, Tree treeOne) throws ImportException { // many trees
        double rootHeight;
        DecimalFormat twoDForm = new DecimalFormat("####0.##");

        Tree[] trees = new Tree[combiTrees];
        double[][] tips = new double[combiTrees][treeOne.getExternalNodeCount() + 1];
        double[] origins = new double[trees.length];

        trees[0] = treeOne;

        for (int t = 1; t < combiTrees; t++) {
            trees[t] = getRandomTree();
            if (trees[t] == null) throw new ImportException("get null random tree");
//            System.out.println(t + " => " + trees[t].toString());
        }

        for (int t = 0; t < trees.length; t++) {
            System.out.println(t + " => " + trees[t]);

            for (int i = 0; i < trees[t].getTaxonCount(); i++) {
                FlexibleNode node = (FlexibleNode) trees[t].getExternalNode(i);
//                System.out.println(node.getTaxon() + " has " + node.getHeight());
                tips[t][Integer.parseInt(node.getTaxon().getId())] = node.getHeight();
            }
            rootHeight = ((FlexibleNode) trees[t].getRoot()).getHeight();
            origins[t] = Double.valueOf(twoDForm.format(rootHeight + 100.0));
            System.out.println("tree " + t + " root height = " + rootHeight + " origin = " + origins[t]);
            System.out.println("\n");
        }


        if (index < 0) {
            printXML(tips[0]);
        } else {
            try {
                outputBDSSXML(curD, index, tips, origins, trees);
//                outputExponetialXML(curD, index, tips, origins, trees);
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        System.out.println("\n");
    }

    private static Tree getRandomTree() {
        Random random = new Random();
        int c = random.nextInt(10) + 1; // [1, 10]
        int tId = random.nextInt(9000) + 1000; // [1000, 10000]

        try {
            System.out.println("randomly get " + tId + "th tree from " + pathInput + c + inputFileName + "\n\n");
            FileReader fileReader = new FileReader(pathInput + c + inputFileName);
            GetDateFromTree getDateFromTree = new GetDateFromTree(fileReader); // many trees
            int totalTrees = 0;

            try {
                while (getDateFromTree.hasTree()) {
                    Tree tree = getDateFromTree.importNextTree();
                    totalTrees++;
                    if (totalTrees == tId) {
                        return tree;
                    }
                }

            } catch (ImportException e) {
                System.err.println("Error Parsing Input Tree: " + e.getMessage());
                return null;
            }

            fileReader.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        return null;
    }

    private static String insetTreeIndex(int treeIndex, String tree) {
        int inserts = 0;
        String newTaxaTree = "";
        StringTokenizer st = new StringTokenizer(tree, ":");
        String tmp;
        while (st.hasMoreTokens()) {
            tmp = st.nextToken();
            if (tmp.endsWith(")") || tmp.endsWith(");")) {
                if (!tmp.endsWith(");")) tmp += ":";
            } else {
                tmp += "t" + treeIndex + ":";
                inserts++;
            }
//            System.out.println (tmp);
            newTaxaTree += tmp;
        }
//        System.out.println ("changed " + inserts + " taxon");
        return newTaxaTree;
    }

    private static void outputBDSSXML(String curD, int index, double[][] tips, double[] origin, Tree[] trees) throws IOException {
        XMLWriter w = writeHeadAndTaxa(curD, index, tips);

        w.flush();
        w.writeText("\t<!-- Stadler et al (2011) : Estimating the basic reproductive number from viral sequence data, Submitted.-->\n" +
                "\t<birthDeathSerialSampling id=\"bdss\" units=\"substitutions\" hasFinalSample=\"false\">\n" +//logTransformed=\"true\">\n" +
                "\t\t<birthRate>\n" +
                "\t\t\t<parameter id=\"bdss.birthRate\" value=\"8.23E-4\" lower=\"0.0\" upper=\"1000.0\"/>\n" +
                "\t\t</birthRate>\n");
        if (isRelativeDeath) {
            w.writeText("\t\t<relativeDeathRate>\n" +
                    "\t\t\t<parameter id=\"bdss.relativeDeathRate\" value=\"0.107\" lower=\"0.0\" upper=\"100.0\"/>\n" +
                    "\t\t</relativeDeathRate>\n");
        } else {
            w.writeText("\t\t<deathRate>\n" +
                    "\t\t\t<parameter id=\"bdss.deathRate\" value=\"9.46e-5\" lower=\"0.0\" upper=\"1000.0\"/>\n" +
                    "\t\t</deathRate>\n");
        }
        w.writeText("\t\t<sampleProbability>\n" +
                "\t\t\t<parameter id=\"bdss.sampleProbability\" value=\"0.01\" lower=\"0.0\" upper=\"1.0\"/>\n" +
                "\t\t</sampleProbability>\n" +
                "\t\t<psi>\n" +
                "\t\t\t<parameter id=\"bdss.psi\" value=\"2.78E-4\" lower=\"0.0\" upper=\"100.0\"/>\n" +
                "\t\t</psi>\n" +
                "\t\t<origin>\n" +
                "\t\t\t<parameter id=\"bdss.origin\" value=\"" + origin[0] + "\" lower=\"0.0\" upper=\"14000.0\"/>\n" +
                "\t\t</origin>\n" +
                "\t\t<sampleBecomesNonInfectiousProb>\n" +
                "\t\t\t<parameter id=\"bdss.r\" value=\"1.0\"/>\n" +
                "\t\t</sampleBecomesNonInfectiousProb>\n" +
                "\t</birthDeathSerialSampling>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t<!-- Stadler et al (2011) : Estimating the basic reproductive number from viral sequence data, Submitted.-->\n" +
                    "\t<birthDeathSerialSampling id=\"bdss" + tree + "\" units=\"substitutions\" hasFinalSample=\"false\">\n" +//logTransformed=\"true\">\n" +
                    "\t\t<birthRate>\n" +
                    "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
                    "\t\t</birthRate>\n");
            if (isRelativeDeath) {
                w.writeText("\t\t<relativeDeathRate>\n" +
                        "\t\t\t<parameter idref=\"bdss.relativeDeathRate\"/>\n" +
                        "\t\t</relativeDeathRate>\n");
            } else {
                w.writeText("\t\t<deathRate>\n" +
                        "\t\t\t<parameter idref=\"bdss.deathRate\"/>\n" +
                        "\t\t</deathRate>\n");
            }
            w.writeText("\t\t<sampleProbability>\n" +
                    "\t\t\t<parameter idref=\"bdss.sampleProbability\"/>\n" +
                    "\t\t</sampleProbability>\n" +
                    "\t\t<psi>\n" +
                    "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
                    "\t\t</psi>\n" +
                    "\t\t<origin>\n" +
                    "\t\t\t<parameter id=\"bdss" + tree + ".origin\" value=\"" + origin[tree - 1] + "\" lower=\"0.0\" upper=\"14000.0\"/>\n" +
                    "\t\t</origin>\n" +
                    "\t\t<sampleBecomesNonInfectiousProb>\n" +
                    "\t\t\t<parameter idref=\"bdss.r\"/>\n" +
                    "\t\t</sampleBecomesNonInfectiousProb>\n" +
                    "\t</birthDeathSerialSampling>\n");
        }
        w.flush();

//        w.writeText("\t<RPNcalculator id=\"mur\">\n" +
//                "\t\t<variable name=\"b\">\n" +
//                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<variable name=\"d\">\n" +
//                "\t\t\t<parameter idref=\"bdss.deathRate\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<expression name=\"mur\">d b /</expression> <!--  d/b  -->\n" +
//                "\t</RPNcalculator>\n" +
//                "\t<RPNcalculator id=\"R0\">\n" +
//                "\t\t<variable name=\"b\">\n" +
//                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<variable name=\"d\">\n" +
//                "\t\t\t<parameter idref=\"bdss.deathRate\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<variable name=\"s\">\n" +
//                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<variable name=\"r\">\n" +
//                "\t\t\t<parameter idref=\"bdss.r\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<expression name=\"R0\">b d s r * + /</expression> <!--  b/(d+s*r) -->\n" +
//                "\t</RPNcalculator>\n");
//        w.writeText("\t<RPNcalculator id=\"td\">\n" +
//                "\t\t<variable name=\"d\">\n" +
//                "\t\t\t<parameter idref=\"bdss.deathRate\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<variable name=\"s\">\n" +
//                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<expression name=\"td\">d s + </expression> <!--  d + s  -->\n" +
//                "\t</RPNcalculator>\n");
//        w.writeText("\t<RPNcalculator id=\"expp\">\n" +
//                "\t\t<variable name=\"p\">\n" +
//                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t</variable>\n" +
//                "\t\t<expression name=\"expp\">\n" +
//                "\t\t\tp exp\n" +
//                "\t\t</expression>\n" +
//                "\t</RPNcalculator>\n");
//        for (int tree = 2; tree <= combiTrees; tree++) {
//            w.writeText("\t<RPNcalculator id=\"R0" + tree + "\">\n" +
//                    "\t\t<variable name=\"b\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".birthRate\"/>\n" +
//                    "\t\t</variable>\n" +
//                    "\t\t<variable name=\"d\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".deathRate\"/>\n" +
//                    "\t\t</variable>\n" +
//                    "\t\t<variable name=\"s\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".psi\"/>\n" +
//                    "\t\t</variable>\n" +
//                    "\t\t<variable name=\"r\">\n" +
//                    "\t\t\t<parameter idref=\"bdss.r\"/>\n" +
//                    "\t\t</variable>\n" +
//                    "\t\t<expression name=\"R0" + tree + "\">\n" +
//                    "\t\t\tb b d * s r * + /\n" +
//                    "\t\t</expression>\n" +
//                    "\t</RPNcalculator>\n");
//        }

        w.flush();
        w.writeText("\n" +
                "\t<!-- Generate a random starting tree under the coalescent process      -->\n" +
                "\t<newick id=\"startingTree\">\n");
        w.write(trees[0].toString());
        w.writeText("\n" + "\t</newick>\n");
        w.flush();
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<!-- Generate a random starting tree under the coalescent process      -->\n" +
                    "\t<newick id=\"startingTree" + tree + "\">\n");
            w.write(insetTreeIndex(tree, trees[tree - 1].toString()));
            w.writeText("\n" + "\t</newick>\n");
            w.flush();
        }

        w.writeText("\n" +
                "\t<!-- Generate a tree model                                                   -->\n" +
                "\t<treeModel id=\"treeModel\">\n" +
                "\t\t<coalescentTree idref=\"startingTree\"/>\n" +
                "\t\t<rootHeight>\n" +
                "\t\t\t<parameter id=\"treeModel.rootHeight\"/>\n" +
                "\t\t</rootHeight>\n" +
                "\t\t<nodeHeights internalNodes=\"true\">\n" +
                "\t\t\t<parameter id=\"treeModel.internalNodeHeights\"/>\n" +
                "\t\t</nodeHeights>\n" +
                "\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                "\t\t\t<parameter id=\"treeModel.allInternalNodeHeights\"/>\n" +
                "\t\t</nodeHeights>\n" +
                "\n" +
                "\t\t<!-- END Tip date sampling                                                   -->\n" +
                "\t</treeModel>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<!-- Generate a tree model                                                   -->\n" +
                    "\t<treeModel id=\"treeModel" + tree + "\">\n" +
                    "\t\t<coalescentTree idref=\"startingTree" + tree + "\"/>\n" +
                    "\t\t<rootHeight>\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".rootHeight\"/>\n" +
                    "\t\t</rootHeight>\n" +
                    "\t\t<nodeHeights internalNodes=\"true\">\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".internalNodeHeights\"/>\n" +
                    "\t\t</nodeHeights>\n" +
                    "\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".allInternalNodeHeights\"/>\n" +
                    "\t\t</nodeHeights>\n" +
                    "\n" +
                    "\t\t<!-- END Tip date sampling                                                   -->\n" +
                    "\t</treeModel>\n");

        }
        w.flush();

        w.writeText("\n" +
                "\t<!-- Generate a speciation likelihood for Yule or Birth Death                -->\n" +
                "\t<speciationLikelihood id=\"speciation\">\n" +
                "\t\t<model>\n" +
                "\t\t\t<birthDeathSerialSampling idref=\"bdss\"/>\n" +
                "\t\t</model>\n" +
                "\t\t<speciesTree>\n" +
                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t</speciesTree>\n" +
                "\t</speciationLikelihood>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<!-- Generate a speciation likelihood for Yule or Birth Death                -->\n" +
                    "\t<speciationLikelihood id=\"speciation" + tree + "\">\n" +
                    "\t\t<model>\n" +
                    "\t\t\t<birthDeathSerialSampling idref=\"bdss" + tree + "\"/>\n" +
                    "\t\t</model>\n" +
                    "\t\t<speciesTree>\n" +
                    "\t\t\t<treeModel idref=\"treeModel" + tree + "\"/>\n" +
                    "\t\t</speciesTree>\n" +
                    "\t</speciationLikelihood>\n");
        }

        w.writeText("\n" + "\t<!-- Define operators                                                        -->\n" +
                "\t<operators id=\"operators\">\n" +
//                "\t\t<randomWalkOperator windowSize=\"1.0\" weight=\"10\">\n" +
//                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
//                "\t\t</randomWalkOperator>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
                "\t\t</scaleOperator>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n");
        if (isRelativeDeath) {
            w.writeText("\t\t\t<parameter idref=\"bdss.relativeDeathRate\"/>\n");
        } else {
            w.writeText("\t\t\t<parameter idref=\"bdss.deathRate\"/>\n");
        }
        w.writeText("\t\t</scaleOperator>\n");
        if (!isRelativeDeath) {
            w.writeText("\t\t<upDownOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                    "\t\t\t<up>\n" +
                    "\t\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
                    "\t\t\t\t<parameter idref=\"bdss.deathRate\"/>\n" +
                    "\t\t\t</up>\n" +
                    "\t\t\t<down/>\n" +
                    "\t\t</upDownOperator>\n");
        }
//                "\t\t<randomWalkOperator windowSize=\"1.0\" weight=\"10\">\n" +
//                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t</randomWalkOperator>\n" +
//                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
//                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t</scaleOperator>\n"+
        w.writeText("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                "\t\t\t<parameter idref=\"bdss.origin\"/>\n" +
                "\t\t</scaleOperator>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
//                    "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".birthRate\"/>\n" +
//                    "\t\t</scaleOperator>\n" +
//                    "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".deathRate\"/>\n" +
//                    "\t\t</scaleOperator>\n" +
//                    "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".psi\"/>\n" +
//                    "\t\t</scaleOperator>\n" +
                    "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                    "\t\t\t<parameter idref=\"bdss" + tree + ".origin\"/>\n" +
                    "\t\t</scaleOperator>\n");
        }
        w.writeText("\n" + "\t</operators>");

        w.flush();
        w.writeText("\n" +
                "\t<!-- Define MCMC                                                             -->\n" +
                "\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\">\n" +
                "\t\t<posterior id=\"posterior\">\n" +
                "\t\t\t<prior id=\"prior\">\n" +
                "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1000.0\">\n" +
                "\t\t\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
                "\t\t\t\t</uniformPrior>\n");
//                "\t\t\t\t<oneOnXPrior>\n" +
//                "\t\t\t\t\t<parameter idref=\"bdss.birthRate\"/>\n" +
//                "\t\t\t\t</oneOnXPrior>\n" +
        if (isRelativeDeath) {
            w.writeText("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"100.0\">\n" +
                    "\t\t\t<parameter idref=\"bdss.relativeDeathRate\"/>\n");
        } else {
            w.writeText("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1000.0\">\n" +
                    "\t\t\t<parameter idref=\"bdss.deathRate\"/>\n");
        }
        w.writeText("\t\t\t\t</uniformPrior>\n"
//                +
//                "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"100.0\">\n" +
//                "\t\t\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t\t\t</uniformPrior>\n"
        );
//                "\t\t\t\t<oneOnXPrior>\n" +
//                "\t\t\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
//                "\t\t\t\t</oneOnXPrior>\n" +
//        for (int tree = 2; tree <= combiTrees; tree++) {
//            w.writeText("\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"100000.0\">\n" +
//                    "\t\t\t\t\t<parameter idref=\"bdss" + tree + ".birthRate\"/>\n" +
//                    "\t\t\t\t</uniformPrior>\n" +
//                    "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1000.0\">\n" +
//                    "\t\t\t\t\t<parameter idref=\"bdss" + tree + ".deathRate\"/>\n" +
//                    "\t\t\t\t</uniformPrior>\n" +
//                    "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"100.0\">\n" +
//                    "\t\t\t\t\t<parameter idref=\"bdss" + tree + ".psi\"/>\n" +
//                    "\t\t\t\t</uniformPrior>\n");
//        }
        w.writeText("\n" +
                "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.7976931348623157E308\">\n" +
                "\t\t\t\t\t<parameter idref=\"bdss.origin\"/>\n" +
                "\t\t\t\t</uniformPrior>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.7976931348623157E308\">\n" +
                    "\t\t\t\t\t<parameter idref=\"bdss" + tree + ".origin\"/>\n" +
                    "\t\t\t\t</uniformPrior>\n");
        }

        w.writeText("\n" +
                "\t\t\t</prior>\n" +
                "\t\t\t<likelihood id=\"likelihood\">\n" +
                "\t\t\t\t<speciationLikelihood idref=\"speciation\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t\t<speciationLikelihood idref=\"speciation" + tree + "\"/>\n");
        }
        w.writeText("\t\t\t</likelihood>\n" +
                "\t\t</posterior>\n" +
                "\t\t<operators idref=\"operators\"/>\n");

        w.flush();
        w.writeText("\n" +
                "\t\t<!-- write log to screen                                                     -->\n" +
                "\t\t<log id=\"screenLog\" logEvery=\"100000\">\n" +
                "\t\t\t<column label=\"Posterior\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"Prior\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<prior idref=\"prior\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"speciation\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<likelihood idref=\"speciation\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"rootHeight\" sf=\"6\" width=\"12\">\n" +
                "\t\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n");
        if (isRelativeDeath) {
            w.writeText("\t\t\t<parameter idref=\"bdss.relativeDeathRate\"/>\n");
        } else {
            w.writeText("\t\t\t<parameter idref=\"bdss.deathRate\"/>\n");
        }
        w.writeText("\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
                "\t\t\t<parameter idref=\"bdss.r\"/>\n" 
//                "\t\t\t<RPNcalculator idref=\"mur\"/>\n" +
//                "\t\t\t<RPNcalculator idref=\"R0\"/>\n" +
//                "\t\t\t<RPNcalculator idref=\"td\"/>\n"
        );
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t<parameter idref=\"treeModel" + tree + ".rootHeight\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".birthRate\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".deathRate\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".psi\"/>\n" +
//                    "\t\t\t<RPNcalculator idref=\"R0" + tree + "\"/>\n" +
                    "\t\t\t<parameter idref=\"bdss" + tree + ".origin\"/>\n");
        }

        w.writeText("\t\t</log>\n" +
                "\n" +
                "\t\t<!-- write log to file                                                       -->\n" +
                "\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"T" + curD + "_" + Integer.toString(index) + ".log\" overwrite=\"false\">\n" +
                "\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t<prior idref=\"prior\"/>\n" +
                "\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t<parameter idref=\"bdss.birthRate\"/>\n");
        if (isRelativeDeath) {
            w.writeText("\t\t\t<parameter idref=\"bdss.relativeDeathRate\"/>\n");
        } else {
            w.writeText("\t\t\t<parameter idref=\"bdss.deathRate\"/>\n");
        }
        w.writeText("\t\t\t<parameter idref=\"bdss.sampleProbability\"/>\n" +
                "\t\t\t<parameter idref=\"bdss.psi\"/>\n" +
                "\t\t\t<parameter idref=\"bdss.origin\"/>\n" +
                "\t\t\t<parameter idref=\"bdss.r\"/>\n"
//                "\t\t\t<RPNcalculator idref=\"mur\"/>\n" +
//                "\t\t\t<RPNcalculator idref=\"R0\"/>\n" +
//                "\t\t\t<RPNcalculator idref=\"td\"/>\n"
        );
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t<parameter idref=\"treeModel" + tree + ".rootHeight\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".birthRate\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".deathRate\"/>\n" +
//                    "\t\t\t<parameter idref=\"bdss" + tree + ".psi\"/>\n" +
//                    "\t\t\t<RPNcalculator idref=\"R0" + tree + "\"/>\n" +
                    "\t\t\t<parameter idref=\"bdss" + tree + ".origin\"/>\n");
        }

        w.writeText("\n" +
                "\t\t\t<speciationLikelihood idref=\"speciation\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t\t<speciationLikelihood idref=\"speciation" + tree + "\"/>\n");
        }

        w.writeText("\t\t</log>\n" +
                "\n" +
//                "\t\t<!-- write tree log to file                                                  -->\n" +
//                "\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"T" + curD + "_" + Integer.toString(index) + ".trees\" sortTranslationTable=\"true\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t\t<posterior idref=\"posterior\"/>\n" +
//                "\t\t</logTree>\n" +
                "\t</mcmc>\n" +
                "\t<report>\n" +
                "\t\t<property name=\"timer\">\n" +
                "\t\t\t<mcmc idref=\"mcmc\"/>\n" +
                "\t\t</property>\n" +
                "\t</report>\n" +
                "</beast>\n");

        w.flush();
        w.close();
    }

    private static void outputExponetialXML(String curD, int index, double[][] tips, double[] origin, Tree[] trees) throws IOException {
        XMLWriter w = writeHeadAndTaxa(curD, index, tips);

        w.flush();
        w.writeText("\n" +
                "\t<exponentialGrowth id=\"exponential\" units=\"years\">\n" +
                "\t\t<populationSize>\n" +
                "\t\t\t<parameter id=\"exponential.popSize\" value=\"100.0\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
                "\t\t</populationSize>\n" +
                "\t\t<growthRate>\n" +
                "\t\t\t<parameter id=\"exponential.growthRate\" value=\"4.50E-4\" lower=\"-Infinity\" upper=\"Infinity\"/>\n" +
                "\t\t</growthRate>\n" +
                "\t</exponentialGrowth>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<exponentialGrowth id=\"exponential" + tree + "\" units=\"years\">\n" +
                    "\t\t<populationSize>\n" +
                    "\t\t\t<parameter id=\"exponential" + tree + ".popSize\" value=\"100.0\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
                    "\t\t</populationSize>\n" +
                    "\t\t<growthRate>\n" +
                    "\t\t\t<parameter idref=\"exponential.growthRate\"/>\n" +
                    "\t\t</growthRate>\n" +
                    "\t</exponentialGrowth>\n");
        }

        w.flush();
        w.writeText("\n" +
                "\t<!-- Generate a random starting tree under the coalescent process      -->\n" +
                "\t<newick id=\"startingTree\">\n");
        w.write(trees[0].toString());
        w.writeText("\n" + "\t</newick>\n");
        w.flush();
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<!-- Generate a random starting tree under the coalescent process      -->\n" +
                    "\t<newick id=\"startingTree" + tree + "\">\n");
            w.write(insetTreeIndex(tree, trees[tree - 1].toString()));
            w.writeText("\n" + "\t</newick>\n");
            w.flush();
        }

        w.writeText("\n" +
                "\t<!-- Generate a tree model                                                   -->\n" +
                "\t<treeModel id=\"treeModel\">\n" +
                "\t\t<coalescentTree idref=\"startingTree\"/>\n" +
                "\t\t<rootHeight>\n" +
                "\t\t\t<parameter id=\"treeModel.rootHeight\"/>\n" +
                "\t\t</rootHeight>\n" +
                "\t\t<nodeHeights internalNodes=\"true\">\n" +
                "\t\t\t<parameter id=\"treeModel.internalNodeHeights\"/>\n" +
                "\t\t</nodeHeights>\n" +
                "\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                "\t\t\t<parameter id=\"treeModel.allInternalNodeHeights\"/>\n" +
                "\t\t</nodeHeights>\n" +
                "\n" +
                "\t\t<!-- END Tip date sampling                                                   -->\n" +
                "\t</treeModel>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<!-- Generate a tree model                                                   -->\n" +
                    "\t<treeModel id=\"treeModel" + tree + "\">\n" +
                    "\t\t<coalescentTree idref=\"startingTree" + tree + "\"/>\n" +
                    "\t\t<rootHeight>\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".rootHeight\"/>\n" +
                    "\t\t</rootHeight>\n" +
                    "\t\t<nodeHeights internalNodes=\"true\">\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".internalNodeHeights\"/>\n" +
                    "\t\t</nodeHeights>\n" +
                    "\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">\n" +
                    "\t\t\t<parameter id=\"treeModel" + tree + ".allInternalNodeHeights\"/>\n" +
                    "\t\t</nodeHeights>\n" +
                    "\n" +
                    "\t\t<!-- END Tip date sampling                                                   -->\n" +
                    "\t</treeModel>\n");

        }
        w.flush();

        w.writeText("\n" +
                "\t<coalescentLikelihood id=\"coalescent\">\n" +
                "\t\t<model>\n" +
                "\t\t\t<exponentialGrowth idref=\"exponential\"/>\n" +
                "\t\t</model>\n" +
                "\t\t<populationTree>\n" +
                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t</populationTree>\n" +
                "\t</coalescentLikelihood>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<coalescentLikelihood id=\"coalescent" + tree + "\">\n" +
                    "\t\t<model>\n" +
                    "\t\t\t<exponentialGrowth idref=\"exponential" + tree + "\"/>\n" +
                    "\t\t</model>\n" +
                    "\t\t<populationTree>\n" +
                    "\t\t\t<treeModel idref=\"treeModel" + tree + "\"/>\n" +
                    "\t\t</populationTree>\n" +
                    "\t</coalescentLikelihood>\n");
        }

        w.writeText("\n" + "\t<!-- Define operators                                                        -->\n" +
                "\t<operators id=\"operators\">\n" +
                "\t<randomWalkOperator windowSize=\"1.0\" weight=\"10\">\n" +
                "\t\t<parameter idref=\"exponential.growthRate\"/>\n" +
                "\t</randomWalkOperator>\n" +
                "\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                "\t\t<parameter idref=\"exponential.popSize\"/>\n" +
                "\t</scaleOperator>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\n" +
                    "\t<scaleOperator scaleFactor=\"0.75\" weight=\"10\">\n" +
                    "\t\t<parameter idref=\"exponential" + tree + ".popSize\"/>\n" +
                    "\t</scaleOperator>\n");
        }
        w.writeText("\n" + "\t</operators>");

        w.flush();
        w.writeText("\n" +
                "\t<!-- Define MCMC                                                             -->\n" +
                "\t<mcmc id=\"mcmc\" chainLength=\"10000000\" autoOptimize=\"true\">\n" +
                "\t\t<posterior id=\"posterior\">\n" +
                "\t\t\t<prior id=\"prior\">\n" +
                "\t\t\t\t<oneOnXPrior>\n" +
                "\t\t\t\t\t<parameter idref=\"exponential.popSize\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t\t\t<parameter idref=\"exponential" + tree + ".popSize\"/>\n");
        }
        w.writeText("\t\t\t\t</oneOnXPrior>\n" +
                "\n" +
                "\t\t\t\t<laplacePrior mean=\"0.0010\" scale=\"2.0467423048835964E-4\">\n" +
                "\t\t\t\t\t<parameter idref=\"exponential.growthRate\"/>\n" +
                "\t\t\t\t</laplacePrior>\n" +
                "\t\t\t\t<coalescentLikelihood idref=\"coalescent\"\n/>");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t\t<coalescentLikelihood idref=\"coalescent" + tree + "\"/>\n");
        }

        w.writeText("\n" +
                "\t\t\t</prior>\n" +
                "\t\t</posterior>\n" +
                "\t\t<operators idref=\"operators\"/>\n");

        w.flush();
        w.writeText("\n" +
                "\t\t<!-- write log to screen                                                     -->\n" +
                "\t\t<log id=\"screenLog\" logEvery=\"100000\">\n" +
                "\t\t\t<column label=\"Posterior\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"Prior\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<prior idref=\"prior\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"coalescentLikelihood\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"rootHeight\" sf=\"6\" width=\"12\">\n" +
                "\t\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<parameter idref=\"exponential.growthRate\"/>\n" +
                "\t\t\t<parameter idref=\"exponential.popSize\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t<parameter idref=\"treeModel" + tree + ".rootHeight\"/>\n" +
                    "\t\t\t<parameter idref=\"exponential" + tree + ".popSize\"/>\n");
        }

        w.writeText("\t\t</log>\n" +
                "\n" +
                "\t\t<!-- write log to file                                                       -->\n" +
                "\t\t<log id=\"fileLog\" logEvery=\"1000\" fileName=\"E" + curD + "_" + Integer.toString(index) + ".log\" overwrite=\"false\">\n" +
                "\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t<prior idref=\"prior\"/>\n" +
                "\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t<parameter idref=\"exponential.growthRate\"/>\n" +
                "\t\t\t<parameter idref=\"exponential.popSize\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t<parameter idref=\"treeModel" + tree + ".rootHeight\"/>\n" +
                    "\t\t\t<parameter idref=\"exponential" + tree + ".popSize\"/>\n");
        }

        w.writeText("\n" +
                "\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t\t\t\t<coalescentLikelihood idref=\"coalescent" + tree + "\"/>\n");
        }

        w.writeText("\t\t</log>\n" +
                "\n" +
//                "\t\t<!-- write tree log to file                                                  -->\n" +
//                "\t\t<logTree id=\"treeFileLog\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"E" + curD + "_" + Integer.toString(index) + ".trees\" sortTranslationTable=\"true\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t\t<posterior idref=\"posterior\"/>\n" +
//                "\t\t</logTree>\n");
//        for (int tree = 2; tree <= combiTrees; tree++) {
//            w.writeText("\t\t<logTree id=\"treeFileLog" + tree + "\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"E" + curD + "_" + Integer.toString(index) + "_" + tree + ".trees\" sortTranslationTable=\"true\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel" + tree + "\"/>\n" +
//                "\t\t\t<posterior idref=\"posterior\"/>\n" +
//                "\t\t</logTree>\n" +
                "\t</mcmc>\n" +
                "\t<report>\n" +
                "\t\t<property name=\"timer\">\n" +
                "\t\t\t<mcmc idref=\"mcmc\"/>\n" +
                "\t\t</property>\n" +
                "\t</report>\n" +
                "</beast>\n");

        w.flush();
        w.close();
    }

    private static XMLWriter writeHeadAndTaxa(String curD, int index, double[][] tips) throws IOException {
        String f = path + curD + "/T_" + Integer.toString(index) + ".xml";
        System.out.println("Creating xml : " + f);
        XMLWriter w = new XMLWriter(new BufferedWriter(new FileWriter(new File(f))));
        w.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>\n" + "\n" +
                "<!-- Generated by BEAUTi v1.7.0 Prerelease r3910                             -->\n" +
                "<!--       by Alexei J. Drummond and Andrew Rambaut                          -->\n" +
                "<!--       Department of Computer Science, University of Auckland and        -->\n" +
                "<!--       Institute of Evolutionary Biology, University of Edinburgh        -->\n" +
                "<!--       http://beast.bio.ed.ac.uk/                                        -->\n" +
                "<beast>\n" +
                "\n" +
                "\t<!-- The list of taxa to be analysed (can also include dates/ages).          -->\n" +
                "\t<!-- ntax=" + tips[0].length + "                                                                -->\n" +
                "\t<taxa id=\"taxa\">\n");

        for (int n = 1; n < tips[0].length; n++) {
            w.writeText("\t<taxon id=\"" + n + "\">\n" +
                    "\t\t<date value=\"" + tips[0][n] + "\" direction=\"backwards\" units=\"years\" />\n" +
                    "\t</taxon>\n");
        }
        w.writeText("\t</taxa>\n");
        for (int tree = 2; tree <= combiTrees; tree++) {
            w.writeText("\t<taxa id=\"taxa" + tree + "\">\n");
            for (int n = 1; n < tips[tree - 1].length; n++) {
                w.writeText("\t<taxon id=\"" + n + "t" + tree + "\">\n" +
                        "\t\t<date value=\"" + tips[tree - 1][n] + "\" direction=\"backwards\" units=\"years\" />\n" +
                        "\t</taxon>\n");
            }
            w.writeText("\t</taxa>\n");
        }
        return w;
    }

    private static void printXML(double[] tips) {
        System.out.println("\t<taxa id=\"taxa\">");
        for (int n = 1; n < tips.length; n++) {
            System.out.println("\t<taxon id=\"" + n + "\">");
            System.out.println("\t\t<date value=\"" + tips[n] + "\" direction=\"backwards\" units=\"years\" />");
            System.out.println("\t</taxon>");
        }
        System.out.println("\t</taxa>");
    }

//    static final String testingTree = "((1:1,2:1):1,(3:0.25,4:0.75):1);";
    //
    //    static final String treeImported = "((3:0.3824976747708124,94:3.401966643328916):0.2131455444432433,((63:0.7035650391093435,(((61:0.2412493040961774,(((69:0.8620006828340152,((48:0.5302937984554452,((66:0.49360327162170314,96:0.7203824406192432):0.26805517248648547,(84:1.6647723439430746,18:1.7899504518955753):8.032342341470766E-4):0.3878430419933654):0.0020541955766213427,(((39:0.547175618789725,41:0.7713476736692841):1.5199644442202225,59:2.300348903048058):0.3545780299502965,(5:2.173657685837823,(7:0.4831518556113801,46:0.04259431830686622):0.9936840610288644):0.3828314786413345):0.06944108321445075):0.2916463175592021):0.8511296243775233,(((13:0.21179557713185626,65:0.9050599247224358):1.2063337719436882,(19:0.961325550099684,(71:0.27227635823039753,((72:0.3623541449244776,34:1.1593527329427098):0.4329707059057515,24:1.750349112100648):0.8548087225791183):0.33826729136137956):0.014057115381393093):0.6343929119615161,99:3.7025837488577316):0.543866697708435):0.1836689943905121,((((((82:1.7276842975724649,91:1.607723421568065):0.1933685581812319,((28:0.23706203140931148,((86:0.5529558273083699,67:0.29663788595939145):0.281373008673687,31:0.8632750147330384):0.5341346272330105):0.7280205172816232,73:0.11489207430514847):0.1042334879532878):0.22023308807550723,97:2.0816724934432385):0.16356342886852548,(25:0.85596627446141,81:0.10864654588741773):0.3265316815053567):0.7052085780791693,((((54:1.2513925633251677,15:1.7216422570576408):1.1052980890372737,79:1.4025557656843413):0.4324573365468609,((60:1.2665130997995677,(23:0.9699990116584232,(58:0.17796506940081697,27:0.3614456732505629):1.3967642598709824):0.4421690408682921):0.7888785083571337,((76:0.1711906024989267,49:1.2710678644023476):0.060131414497704094,64:1.2593916884637246):1.5454624508194899):0.24056330854342534):0.015283090096886554,((32:1.3278266797752578,50:1.2629774321495602):0.7777520083215703,88:0.02772431983654222):1.1101361895232151):0.5963623233265598):0.3665009131043293,((((16:1.9309003970287875,(((90:0.5445216783903744,6:0.8548831149642158):0.2741249294450656,30:0.008362111012826023):1.022234909613422,(40:0.7366132749777619,83:0.5070914524648216):1.3875613903422281):0.5001700205337869):0.7233479987224163,((1:0.6003333838940725,98:0.3415505496586113):2.588962722413393,(4:0.2759376488056353,29:0.47002294631420005):0.898552467271085):0.2362157981780597):0.11672800780799664,((10:1.1355322499789322,43:0.5843151524612831):0.701686409717148,52:1.9374978201794926):1.417945520858075):0.7586500833008447,9:0.3178627566272416):0.03155457313484966):0.14868731051955386):0.3667216608428072):0.02390814464028157,(((((35:1.086623144974067,74:0.692186651442158):0.10673574135114072,(78:0.1638180965937629,89:0.250663800794962):0.08265189041229482):1.4629757070659404,93:1.6063646891197325):1.4550877544464176,((37:0.9226166936844796,21:1.044343709701885):1.5676034487217074,100:1.3447223957080423):1.088773395481709):0.08284086233381505,((((56:2.355869340652318,45:0.4032979506669325):0.3330320301280212,((2:0.5357734051244822,53:1.0622382154097905):1.0714020295534217,(57:1.9553678171881557,92:1.9351676011453818):0.3176309095081522):0.33889485298905964):1.0910719751499318,51:0.06358005997870864):0.24119223678720436,(12:1.0464825190142348,8:1.0984225849334481):1.6150168196252994):0.13538891424735233):0.654056534527701):0.68828126623945,((((95:1.5084591517564723,87:1.0618807586778019):0.13101510403314265,75:1.1323197425395333):2.2878458817626433,36:2.864052872624789):0.6509399099742286,(((17:0.06524692359560014,((((33:0.017056678783653467,(77:0.768977532338945,22:0.2066823108679694):0.3345642602125296):0.600193544474839,47:0.6997906421703952):0.5100738160877549,(70:0.5677418730833756,42:1.9551556564280363):1.6169469279243813):0.3307075738681533,20:2.492040542806981):0.525443622448492):0.09433286051003886,((11:0.3113555135611077,(38:0.44389504903251376,44:1.2708469141489145):0.2622045045057788):1.2732398364643553,(26:1.6716254866488542,(55:0.9828398772396151,68:0.21551121655825067):1.106998506688205):0.5526554864260289):1.671928572660197):0.12614845672773534,(80:0.34195834822687576,(85:0.09264904116899952,62:1.713709090294202):0.18592099291105324):0.4792020544951763):0.2727014597498254):0.2650194727743962):2.4225764992771586):0.12072071192465827,14:0.47972787703636754):2.757605090929111);";
    static final boolean isRelativeDeath = true;
    static final int combiTrees = 1;
    // T2 sampleBecomesNonInfectiousProb = 1; T1 sampleBecomesNonInfectiousProb = 0
    static final String path = "C:/Users/dxie004/Documents/HIV Swiss Cohort study/Tanjia Birth-Death/20110718/brd1/";
    //    static final String path = "/Users/dxie004/BirthDeath/20110616/exp10/";
    //    static final String curD = "10";
    static final String pathInput = "C:/Users/dxie004/Documents/HIV Swiss Cohort study/Tanjia Birth-Death/20110718/T/";
    //    static final String pathInput = "/Users/dxie004/BirthDeath/20110616/T/";
    static final String inputFileName = "/T100.trees";
//    static final String[] trees = {
//            /* 1 */  "(((100:0.2759295941079598,((6:0.1780547126475307,90:0.2773585176913773):0.4324798358805104,(77:0.4366037591568581,(((60:0.06114593641368682,(((45:1.1858532193532096,69:0.4857462933631136):1.0703793675170519,14:1.0080965286106625):0.7531830236235892,(((72:0.4734480373733323,13:1.6414737257868168):0.31345510548895206,9:1.2799154503620205):0.36580160160488395,55:0.6521307508046301):0.13222292659177048):0.6843395008777939):0.08358823875237409,(63:0.1173962764017622,(((((83:0.06388497493919498,(87:0.06971104871286249,96:0.11513746603637598):0.2617677248167065):0.23523100742624647,11:2.5488642344318637):0.5709241701969603,(((59:0.13001207544150828,4:0.534509372507409):0.9065013929723875,(((78:0.4300613180395787,(54:0.5881259163967514,37:0.09161743362459873):0.10007317973276297):0.3581939195808168,71:0.8417847156069083):0.6760435746540918,56:0.34728849185761557):0.0011569801790649592):1.6147112046099765,(98:1.8505753493501458,17:0.9280241884545952):1.5062115835297483):0.6173650351656441):0.0175798939167926,40:1.4327849159224866):0.6484216994792886,42:0.5649331356339093):0.044490735903595535):0.10065007060014697):0.7480152893310459,((19:1.0453148233101777,65:0.03982690649499432):0.015852971398985005,((((((((24:0.6426142775167778,49:0.26507022798260027):0.7090057146513142,((41:0.14186057102002847,67:0.026961796503665925):0.36991918712876304,79:0.5815959657493615):1.024603334843216):0.024183855106581387,86:0.08436576474963875):0.3946779683409196,((61:0.21049623087835378,(16:0.22329709973960155,93:0.0980200763866037):0.11081431464717495):0.8467610850802538,31:1.2165550002603938):0.06206931566212637):0.7024059789204573,((18:0.5110908536902983,34:0.21183915417614707):0.882488365662339,73:0.23324499594150505):2.0871456603448966):0.44201612676651125,99:0.455825827781966):0.638245486290975,(39:0.07658437481359659,58:0.1764000007473152):0.8470799448815418):0.13252970716211898,53:0.5440963965993131):0.7116321144589479):0.04446103714470162):0.5570925454659523):0.8349620719218072):0.39537618606560354):0.23246534833540267,((70:0.1996975744937446,((30:0.25043743640342697,((((89:1.6127193503514903,(82:0.14048148155413376,57:1.0784083186686084):0.38780755293423397):0.24644757705303455,80:0.24722619589322203):0.06489043744710532,7:1.6582683265776672):0.4481749866378619,(2:1.8804071294900604,3:0.023035004075473875):0.7446717171494495):0.2814795274087576):1.0069357276627375,((((((50:0.009904531387953774,(51:0.04974396264227465,((28:0.6164291583816088,33:0.7028750326597675):0.020632702991932206,52:0.2620934650458022):1.6169083530976258):0.11188677321436113):0.5897740532249305,((((84:0.00915957993005767,75:1.403360151297104):0.17953794102191267,81:0.0627513917657474):0.37241766310854807,(((94:0.7136713627250614,10:0.008279153277477747):0.11337143032465669,(68:0.9476218248680841,8:0.46104904327810214):0.3082656519911391):1.0277836373882483,35:1.139900635404104):0.3649293864599641):0.04053080400729936,62:0.43372394838795225):0.8450036334702116):0.21676992831331532,97:0.3937210281260244):0.6217923509933945,66:0.3065474310864129):0.8170463207263605,44:0.026894004345598965):0.015539344835022284,((20:0.1252746831931999,12:0.2706972912118957):0.23084681536101925,(((21:0.04108670219976496,((26:0.21683708421460812,(32:0.14412489943602447,(((74:0.18560721020572668,(95:0.362902776402263,((25:0.42464905250257207,1:0.5164820341492313):0.5925646746092141,85:0.13636374873279244):0.0031874038009704897):0.0044990893278598065):0.10197883895896687,(43:0.05393610804603166,46:0.926216914294938):0.1084119847356193):0.34337795022535045,29:0.5918636480714469):1.0448119892689047):0.31409259065885076):0.6919113759361757,91:0.7479793351634823):0.6136793015327289):0.27258213498288786,((((((48:0.1306215018764466,(88:0.8711354400164752,36:0.4491587803517647):0.055997269374821856):0.06583165565941829,27:1.0506149946403638):2.364007453784181,5:0.45338514679037667):0.053644616923241184,92:0.7051063242279012):0.35121217569111796,47:0.5874377988158308):0.3714369733787217,(64:0.032098512520678746,15:0.07581713088827158):0.03448907859428907):0.2810624510789914):0.30792407694989077,76:1.1156644058672276):0.39531170115012415):0.08177480035578455):0.5568135150951967):0.14934174255400912):1.1303652342548354,(38:0.2423418362203531,23:0.3173102378974546):0.32247555533508354):0.5340110001664744):0.3446866710113099,22:0.12145978249112854);",
//            /* 2 */  "(41:0.5917472044147383,((((22:0.004822363771718408,12:0.036453899723426275):0.038192624923388596,85:0.005582055959308274):1.6405071172099777,32:0.6601791552306269):0.9614057600012043,(96:0.3000557138093116,(((71:1.4169394049035824,((((91:0.1370559611415103,21:1.0666517759455028):0.3013823716294357,(16:0.0573837354307003,65:0.01605939620417085):0.7395436085286314):0.39351458733091915,((4:0.5766757437033982,49:0.5924370722126595):0.6830223876929804,(63:0.44544320701755735,((78:0.5356237019093407,69:0.39872973970167014):0.22658400724539007,2:0.007295231526302626):0.2619316314185782):0.12247535755605843):0.4613297883047476):0.22135415594459662,((66:0.0909903756180076,79:0.23944532548434083):0.9326451214609349,77:0.08611690710177422):0.4928995981351125):0.35086338774172576):0.42888176198105565,25:0.22372954562462066):1.6561895153049986,((8:0.7252985413015516,(((((19:0.5731481928179369,38:0.5628148179977828):0.7280619924310128,(((51:0.774714515779841,(20:0.30259361812520247,(((45:0.22547818679459386,(15:0.11155583153541332,89:0.0019180047209065854):1.2608754603900336):0.20066176408596847,5:0.011746941347904327):0.1258811004867788,99:0.005171249553354684):0.2856988544773096):0.5919984583930877):0.6680436985029119,62:0.4929500980465096):0.3407459001243027,(56:0.16112858782823203,((((86:0.9999158037021223,(((((6:0.3124148501489347,9:0.39002673818040207):0.15232159499193254,42:0.025343389649686188):0.15783625157397607,(52:0.08252282465787572,1:0.1903249012293786):0.8355165755561378):0.8456377495528147,(29:0.09306183195060469,(36:0.5121393333926862,61:0.6151895034113255):0.13125039698954954):0.36730220486517884):0.22509374921401326,((70:0.548251196265543,(54:1.068022441103647,98:0.3918098219375371):0.17015240619493):0.13477388544955327,76:1.32875931188856):0.5573870865863535):0.23948185129769595):0.660980035309358,((50:0.3183515122501228,43:0.4330103644785115):0.7161229906550692,((46:0.6982485369753215,28:0.35914233093703896):0.9456063575498628,((60:0.15152705573659953,88:0.9471691368291206):0.006131445480609177,57:1.2055930837668711):0.1168845181106466):0.2821811539741428):0.28658450213926434):0.12231958355958028,100:0.16886198799975904):0.39432418425509486,(((35:0.5476378565836209,(((18:0.225425185335262,33:0.08978391582746603):0.006197436593837091,(84:0.23753266096590897,47:0.054635169404989026):0.22815199471594472):0.040090186625441504,((68:0.5202211551347012,72:0.15964331441430202):0.411489869130763,(55:0.20859691653340207,73:0.13553593601346248):0.13666515596971207):0.006818326385002571):0.23567377464165418):0.4714343657989333,48:1.5927409626470213):0.3766517469902353,(67:0.12673772400627192,7:1.1462220199153317):0.974055300069147):1.1623642735878006):1.5495037998203287E-4):0.18315690644391847):0.04635664312022536):0.07334039689344474,59:0.22350091373092962):0.03682424644454185,(((80:0.05236949965089277,94:0.23241336373873311):0.9433465766797506,30:0.037619487047705125):0.14754844432500747,(((64:0.27544712472171096,17:0.32206323453342955):0.10492373292922874,(3:0.04182901501472891,((34:1.2737727909122611,11:0.04061530717404138):0.6145476708779154,58:0.4772657820579229):0.5726791297976788):0.5339313711611564):0.1932764110884433,(53:0.2898114246653787,((23:0.45593536930150247,82:0.8656995404468664):0.3872574897975163,((((90:0.38032240461584155,13:0.47798888209289925):0.17639524849569743,27:0.197037289443029):0.16238512947848016,74:0.026071102004512126):1.3140863224468318,((24:2.258874030528173,92:1.6095849857326585):0.029794673691030393,14:0.502284744602526):0.32014871387072485):0.4332773701450061):0.07368159395796603):0.4255757485182081):0.14154533689937532):0.15314283412876462):0.20316950598011063,((81:0.3780202281433547,44:0.010725176953772397):0.09495371670612895,(40:0.12414949153226029,31:0.5324837044517716):0.3872132895991154):0.5652026433796182):0.13878612969561122):0.5163154959983132,(((83:0.010850489430678234,10:0.039730311543561214):0.00451766365623385,26:0.019161324722043016):1.8474257478546705,((93:0.3326316546893904,75:1.630889359942298):1.3819910060920821,((39:0.49621150559368643,(37:0.44370634124779795,97:0.08068959368103207):0.2850705764937068):0.4555503144811305,(87:0.426393053197235,95:0.6627568262239785):0.5328647905832615):2.6497694981512687):0.41591761250949366):0.34497135700431514):0.008789846631671594):0.9178353805692847):0.6415881181349361):0.4310608834467935);",
//            /* 3 */  "(((((38:0.20945893103632507,67:0.2545479255397307):0.2222080071554533,(((65:0.31774859025521174,50:0.07620494738921835):0.511670460174015,85:0.03794764374746862):0.13373395694374146,60:0.2915010611112834):0.19292070159626284):0.7510849391983871,(42:0.42855486719470903,(((63:0.4582823839850825,(((27:0.06037136613195582,((86:0.8957230982196371,((((47:1.0876476672004083,79:0.7621070044901401):0.3376276785627432,16:0.13878068280908806):0.20702134127963667,((((((64:0.018276586754216506,(((24:0.5740147182619897,51:0.42639104987822896):0.32421320524260455,(44:0.8803456367835127,((57:1.0683310288125438,(58:0.6862531509910688,(8:0.430013497883464,37:0.14125995581964101):0.599229054012099):0.030301804680573197):0.06292735899974655,22:0.10081717347691099):0.10011244948457043):0.09604306051323697):1.3542410309549537,54:1.133612494118391):0.0030180869682210343):1.108182455849681,(73:0.2511326714839113,59:0.7513681943684358):0.6440442650082128):0.008753391773896002,((((4:0.1407154090142544,87:0.10374014119431707):0.1975453197206185,20:0.3765702595001357):0.5738151926793509,((68:0.975351222203803,(12:0.35922406519377004,(43:0.17987010921679436,26:0.3936334594191574):0.2144702754592993):0.2126942646121558):0.284290661461547,(97:0.670930362000043,29:0.45435268212249835):0.647479266319371):0.2948457500708177):0.3701633192136722,19:0.20758256418923215):1.796434872270015):0.03523791698211065,((55:0.22065938446675304,31:0.3281976310359913):1.1731485794688266,(21:0.9203936215475217,(18:0.9405704959476324,(((89:1.3016359089711877,(46:1.1231719777003286,1:1.2997693575898497):1.2145534332520014):0.1414433717737853,(96:0.9467745613694305,13:0.009091421818136602):0.004117634916118718):0.9480455844029891,80:0.21065355514420192):0.12983898260500748):0.08641438208265573):0.03403368275555607):0.02696801550694916):0.9188605033069219,((49:1.8751475635664305,91:0.23997884222228505):0.971583078731602,(72:1.1079609618600315,2:0.6690219812752303):0.7536702933599404):0.7335134063491546):0.9722705746091149,71:0.1998659027383649):0.04282118073808494):0.06924540377044774,(70:0.07782821709723997,((23:0.053955970385401564,(((94:0.1259661733564399,((((7:0.4129984244303433,((53:0.3346791590896182,76:0.19703779443461183):1.3006982251854122,(40:0.3199236700753967,5:0.29240810767027714):0.7437410708062011):0.22459654760464365):0.11259202385590905,17:0.11920963478805868):1.138680469035327,56:1.0923247967265817):0.17109242963366134,((34:0.1590251305241952,98:0.5495814898246647):0.4847274065698195,(77:0.42792955513817676,(48:0.8982496932991122,69:0.028245620612876943):0.1953893988437101):0.7234967683583065):0.7995854614340314):0.17078626357896587):0.031742310394827555,28:0.49961022133586086):0.11585793955393253,(61:0.19723668031090558,81:0.22390092119131655):0.16398971737028267):0.5950136595790352):0.2426565942524741,(74:0.39080856473024106,92:0.09304519111921206):0.020287748400023453):0.4905334846330547):0.7895582476861076):0.885610801674324):0.47533228457303167,(10:0.23847259111848018,14:0.3544437659446791):0.26365579922420856):0.19957937389126634):0.8259354324509038,93:1.1113311393599945):1.0403049312127894,100:0.71418665896379):0.18737578322699733):0.9311386748992216,(41:0.30513975257813186,9:0.1369719748629894):0.6601136387806346):0.6190646695272104,15:0.1294129241627786):0.18830155513736813):0.38891731948782216):0.19794438020097438,45:0.5862610066122542):1.0718124455450653,(((36:0.2913177399704825,((90:0.5552880241882718,((32:0.7479380663619146,95:0.47959514427707184):5.233199806653488E-4,78:1.189268612429517):0.10919740376473719):0.7443342957766994,((25:0.14932169558096042,88:0.5232749980187332):0.40645769335340987,(6:0.5073672639081792,(11:0.5186461551353752,((((83:0.03715113773255563,82:0.2614250212766853):2.175681337194349,(52:0.09423981546368676,99:0.2647423126132731):0.6392657973520794):0.22033697732996993,35:0.3131657223274318):0.7309269197135797,((62:0.45652484324192866,(75:0.029001350347898303,30:0.03606563391990392):0.09356574068423384):0.6597114808400275,39:0.33057992125784175):0.491794075591081):1.8882986408869158):0.5962601572361033):1.9156514806670417):1.901609813062942):0.6254831082284333):1.0822951207746563,(3:0.9429879121409019,(33:1.5282283228151101,66:0.37795674680008773):1.5778275620696398):0.1932596604934531):0.7484342011537404,84:1.3028298578034718):0.07192939096426088);",
//            /* 4 */  "(((63:0.5971341772682925,(((32:0.33088442935805706,42:0.9818064820555488):1.2852789040403396,4:0.13844816725331643):0.9050067937827864,((72:0.11417102495455467,((45:1.093769638861756,(9:0.37013162632819674,73:0.6575430958348033):0.4555728093701473):1.4412695815738972,(92:0.4801417763123417,2:0.5887611487368416):0.22590464841751157):0.2646555154818171):0.41874643609900675,(6:0.5010735799172772,88:0.2158783182343491):0.6600982215947955):0.46104736197928986):0.13061301027676242):0.8481434396452228,(((((((((22:0.6522107792498546,93:0.7995169966454267):0.4113143508752384,87:0.6694193558953427):0.5023266366105756,((13:0.31957985443240466,(61:0.09579431210066092,74:0.41647078615105204):0.3152390654347944):0.548201458839291,30:0.05629404375512892):0.23024850383335083):0.21685794963441918,(((3:0.8419301810112902,64:0.42307005820123145):0.039576328944843864,(53:0.21472032439718214,43:0.516016144807318):0.5522246100593633):0.7936702784986178,35:0.21690341488244647):0.05708285546235148):0.6910012684088946,65:0.04950621815533074):0.6568390165260207,((((27:0.8977893879571093,(36:1.3150826865822502,97:0.2909338557526886):0.17749747011945027):0.35303080406711485,85:0.0762044554733925):0.056372912249887275,95:0.9593645074528676):1.1077705071412955,(33:0.25300787069890873,(25:0.328817766496198,(29:0.04569368473045543,(20:0.6500779516561956,58:0.9473702542902976):0.16562370827308825):0.44200373669096393):0.5891576662646738):0.4490136088204424):0.004221652274960519):0.16607731991498742,12:0.7120751581749354):0.06529281018453803,(70:1.029108169555296,((15:0.013060730855017244,(71:0.14523976087067103,75:1.3879194852363435):0.3333589885555921):0.051095275716666944,((21:0.23379576040460553,5:0.08450396831842438):1.2724933282830246,55:0.14982543268533566):1.5819959456094743):0.04238501767700997):0.42445024487807315):0.9335752130467556,(((((59:0.5583765773710039,82:1.0229814889739464):2.0262284935218986,((((56:0.12181364928154137,((60:0.6216053508855168,80:0.7729964858383543):0.020798946837480514,50:0.8822032898467224):0.435092024954284):0.5210573955778481,68:0.22987542095442626):0.37082246837350796,(11:0.2537773523886817,((66:0.1135105187045009,47:1.4477336353709616):0.08562725773698254,((62:0.7970273485498136,(7:0.11598325830248846,1:0.5036319748241439):1.0526505463359794):0.36590489402030246,(52:0.13340072272101278,91:0.1370546630789482):0.46023881312192794):0.02916959442845113):0.38673537889656706):0.39395868912414267):0.7662491502961601,((37:1.110838344710146,(((38:0.9546481672843476,((69:0.19930238000643707,19:0.6005407545168806):1.1738316084346145,10:0.2551589871022202):0.06215257026410326):0.14250195516227704,(40:0.06886909144807296,((((((17:0.4318672095227283,(28:0.21978697257588475,86:0.03310590083384446):0.5183094602583803):0.023033822015270378,(31:0.17869110898513352,67:0.05718399747186875):0.5752856923994089):0.23278629636613868,(81:0.3724158200364376,79:0.7202517037692553):0.33944175123452003):0.13600438064238518,98:1.3005332208675515):0.08332934278262116,89:0.7225868497783603):0.1119773196570013,(96:1.2696299907945414,((48:0.1594873895676261,((51:0.9125092659214435,34:1.1038362280470857):0.1547607142716534,46:0.020852996045211736):0.029486834734979483):0.007497717851954633,((77:0.452520031779162,54:0.3031274472727311):0.261605600170702,24:0.04227110030653991):0.460488804375494):0.11580824023953551):0.3019265040090098):0.4053199474113922):0.536672565288562):0.011357139898990454,((83:0.9616030616846336,94:0.7024704034658114):0.17522950578763918,99:0.44779756653857117):0.017264833109371658):0.08434130999387612):0.025630671144057704,49:0.19121210472840522):0.45648465294874896):0.7919343107310728):0.10338771978532169,8:0.34674292192684764):0.08267662759021555,(76:0.30439306935274946,((44:0.46079144735569333,16:0.3321601586351881):0.21713822199622346,((((14:1.3256098632440168,100:0.25154989345941625):1.0412638747346932,41:0.5012105641168794):0.5387701436727275,78:0.12391763068654527):0.2886549771279019,90:0.3478513544360764):0.2535068195701742):0.09103915532813911):0.33288272351850434):0.1107114076798883,18:0.33107611762760136):0.0925520120575829):0.11415229255684256):0.8351976614004988,(23:1.3172927167590816,(57:0.16887040433248224,((26:0.5563912527182566,84:1.1354131008126298):2.9402470919902255,39:1.298614026879397):0.13619115668976356):0.8768083947557272):0.1410510185120719);",
//            /* 5 */  "(((((14:0.06909967074591528,18:0.008348526065371509):0.6199385534414614,(((((86:1.3981866513993146,(22:0.4881370076469834,(1:0.4025199081725816,31:0.36911750360764556):0.1667172441271389):0.8618265297835835):0.42764098995838706,82:0.284276599882209):1.0262901396467705,(90:2.1896351150424467,(73:0.05785032663113321,((53:0.20728929479767944,49:0.10029765055438067):0.21544168455893864,27:0.05036465068452878):0.14548392795725418):0.8132561453660818):0.6733123617557313):0.15503141878426785,((((57:0.4338727446686934,33:0.12777429823861874):0.01848542631124328,(45:0.2587084986996797,((20:0.2969951121192853,85:0.36345684248360033):0.11821434720863178,15:0.9610866570072945):0.071474268687171):0.7796089161360698):0.06787867335300124,(60:0.8458907738040824,58:0.05684723683962445):0.5138731178983449):1.0504064182935218,(100:1.3449632993117353,((59:0.21126555094416744,91:0.692367982396022):0.19655888203227412,(74:0.33956540090437537,(67:0.21345619080552436,55:0.8998759367337437):0.5717035145499962):0.28597367669129947):0.3458544418358547):0.4788892574841346):0.10137981987700817):0.49883571366511426,(51:2.273775714547024,13:0.07272411702726611):1.2209874472349247):0.521076756165876):0.003694570964231758,((69:0.04206467181822804,(((((44:0.29905042820576977,80:0.3390417228377617):1.0144317679586936,83:0.050106714860378876):1.1812672337421168,52:1.2531971000364592):0.02358754733485391,70:0.29357585267344):0.27367876652559175,(16:1.5708494969306326,21:1.5455056524561397):0.7831599355851946):0.9764375206539095):0.018156706152054003,((97:1.5148456577494125,71:0.5246690125380193):0.07097065774823141,24:0.00796403982324101):0.020251856131282775):0.09815278519497106):1.7529274318141512,((46:0.4109712659865652,(93:0.570274206150633,((35:0.1961383428478949,17:0.9921071840272186):1.299181504487429,(((((((2:0.6217027083854507,11:0.1891019933086091):0.12856372761375512,94:0.34779730454202784):1.9077807293710576,29:0.17934268531349673):8.384961698224558E-4,((88:0.14024533879559176,((32:1.051013208159703,(((9:0.1058348957743559,41:0.5567579155764354):0.2631883617620079,28:0.6930870352289574):0.9920113293622722,(((4:0.49972450122803524,63:0.18557769376959699):1.6666481051108573E-4,54:0.7661490886443401):0.5409339612231581,42:0.2031473083262434):0.4825296181352956):0.31233225100781103):0.20564436461532498,76:0.07842107233367068):0.014030709994730106):0.21011313986110558,(((19:0.9979308951587191,(10:0.4256739992591998,62:0.6261078135311389):0.0779527165717655):0.1286844933126159,38:0.5942268614558595):0.3731957884853887,79:1.1101712959586396):0.9063905116202473):0.11750975067738212):0.13776239700878667,96:1.2462085602022133):0.12169701169988567,((50:0.19011562885314381,(25:0.11605092714413656,48:0.9037564933255857):0.15582813661458061):0.4920972897237488,(72:0.915622847439193,81:0.8035378655879422):0.03659526697898219):0.2518329086891824):1.2109815932384644,78:0.3592170406682089):0.006774675898328475):0.7605547591150037):0.10538157033163742):0.4057033426217096,((61:0.22777104834421102,((((((98:0.37619755379267783,(43:1.1482551151698017,(66:0.7523278709453546,23:0.6979851883565462):0.8140482163710253):0.6366417305228953):0.9501229297463536,((7:0.3789273599150307,(92:0.2045643510327264,((3:0.3838100072951548,64:0.07928612864369611):0.7354763043701136,77:0.2907372366124954):0.232659104032793):0.007437910192638064):0.6819571188651095,36:0.2964825578677073):1.094749969768622):0.5236485031240594,56:0.10123767910634518):0.5655443898160564,(((40:0.7329209921223252,((((75:0.8640066366663972,30:0.8497918872125507):0.14789516915738643,(47:0.1382885341897122,8:0.2502450578723663):0.1949062043386478):0.48052499220805256,65:0.3523270072322342):0.5731975322542326,(((39:1.2012755402628539,84:0.07172284256342576):0.09610494695209115,(87:0.2257137090959176,89:0.09331172687770861):0.388880126465341):0.06393683028345443,99:0.23184422695864448):0.7365986599714165):0.7341531748128616):0.35563295423312846,34:0.3183948196370685):0.21484112721304394,12:0.5263288233966033):0.2878371746600479):0.5873181279956823,(68:0.8053867295135282,95:0.03773728634007956):0.01777980553155878):0.37222209704579257,(26:0.25226043555482036,5:0.0863376023858109):2.2940877664535617):0.04529443408887257):0.023232231434962003,6:0.1161858917291072):0.06709043731940145):0.36940580849439364):0.496188641961111,37:0.4496901978141734);",
//            /* 6 */  "((((66:1.8051857271854743,(27:0.6631995440727945,87:0.9493455897736638):0.14058697866630343):2.044339644575941,(91:0.24714044094438936,(52:0.010082026572115765,(46:0.3384107979262563,((((69:0.2563081951510289,(56:0.19209850379167603,(85:1.1487014739801298,44:0.028141966769021032):0.07278053255357575):0.8060945140472255):1.7031111101108491,((39:0.8121885192358218,((93:0.3212728375101488,((62:1.0714318480929954,60:1.191316504799483):0.2899545086962476,((6:0.05247359857018419,70:0.08016642111032357):1.6766912016292168,(((16:0.2710703447059435,33:0.8031429205917207):0.1995980164676372,(((18:0.3137521647733166,(51:0.17919542687575435,81:0.532183293809031):0.2211625695998567):0.246579225138021,75:0.3036524716677096):0.033318344448690196,1:1.1498681036277671):0.07494081500081617):0.9030096435285131,((78:0.2623254129121191,(((10:0.44920236229616384,43:0.7348634314748519):0.07873907664155344,96:0.34611492300694446):0.025153277777600014,2:0.15288135450516482):0.2297345940020432):0.524541385566262,97:0.2718940196532993):0.009099201418466762):0.35392756731770936):0.451259495223975):0.30074234384011334):0.30932973391370044,94:0.5441203148574121):0.08944899410833962):0.33725457603508646,(98:0.24160091351265178,28:0.5224635765351042):0.02100141154922408):0.01994495313557687):0.14922615986750776,(((((((36:0.8640770396416508,34:0.8658995356131345):0.6662363133067366,(((((3:0.004470690507340169,((86:0.11696262326765483,8:0.01946380269361858):0.6409740060528389,92:0.0571161024376301):0.46434119972766363):0.2617847150605239,(26:0.4659738621260723,(53:0.465643284428392,(7:0.6556440821438937,50:0.8220993308814497):0.7847585965832209):0.2315420749447039):0.025413585150737328):0.4495480716404432,(42:0.5066019935710511,(24:0.16441848539873627,63:0.3109270929376078):0.6560804110885441):0.014610292868232122):0.3575956742976465,(31:0.1824878019912024,((((22:1.4788525055241608,29:0.10692888087324404):0.019846907901590827,100:0.07709417898746329):0.39904728032626813,48:0.20800318210979185):0.20974827021803932,40:0.16354565951487965):0.07338433934596988):0.16532098465160328):0.10854076817052816,(67:0.14973524550154016,(17:2.004017017661622,15:1.2896217656895907):0.08542696167592956):0.39569362961755905):0.18067539030754753):0.28141398648723825,(72:0.8608299131701789,(12:0.6125151396831114,(4:0.7444430438921765,(89:1.1991371200512142,(83:0.9377077160635627,64:0.6926152090937752):0.4659103333472263):0.010920711219980284):0.25725112025023833):0.5928759990406482):0.26538851971680355):0.18058095987412592,(((41:0.025317421500515414,(38:0.5177788914684798,74:0.3597803759840962):0.6467550787164872):0.6881377147089338,73:0.04319216737762743):0.3101215161205584,((54:0.1747510148532217,(99:0.051757307412777176,21:0.12022292279084779):0.040674674657433296):1.2419012411135124,(13:0.1889171299342245,(47:0.09367603200759,65:0.21405571371391408):0.01476240978756449):0.03584046369761684):0.5191383001939949):1.2315875122814122):0.043157274310656746,30:0.5900824386307497):0.3023574350702356,59:0.2564913261170587):0.3152391911170924,(((90:0.42740753167627976,((32:0.15213136211575362,20:0.8796148450215902):0.005020649011782563,61:0.5189674781385305):0.08194851288547556):0.366599449923938,49:1.5314569288514526):0.6720365004224296,((((88:0.5222590371125548,82:0.134764595276595):0.01664471612099061,19:1.3790498661692407):0.19192387207098038,(68:0.045178879969167696,(((84:0.09408886605284161,45:0.20025711686633096):0.04932228749019729,(9:0.22332416564681823,25:0.3869182826178045):0.28086443384539994):0.5084915763967757,58:0.08751403211001607):0.31424869150308843):0.03523442710528646):0.8560567696526316,11:0.8759233768998129):0.4475166845973648):0.7582328044867799):0.045436345892345464):0.27220466278639677,(14:0.08724136422846618,((5:0.9911226622901959,76:0.7920881888799691):0.3244238195974183,37:0.571334180212788):0.19443954198784086):0.572737330471933):0.21025852084237417):0.0138933027769097):1.3000766438965536):0.9822840947855207):1.121349654877636,(35:0.2639216882487281,(57:0.12944428453131707,55:0.37884342879439004):2.528599883295273):0.8653933012996609):0.44010115154763696,(77:0.3110662001971898,(((95:0.6718413317050178,80:0.1938275652432324):1.676089846004932,23:0.40600756750884504):0.39793949171176024,(71:0.2766950876509142,79:0.021716057566452562):0.755124877666371):0.0559359693561543):0.5367300018452985);",
//            /* 7 */  "(((84:0.5218034996774268,(34:0.7248627225527633,58:1.106941990535364):0.7769048713276465):0.15318323458264604,(((78:0.05152105841888144,47:0.1494774579665945):0.22485003835899597,(((74:0.6780468472435803,49:0.08067535008314763):0.2969275939747762,(((98:0.872811312943548,(((((7:0.5189876728272584,(22:0.8091373750894277,(56:0.012720215475482965,95:0.25553876262503517):0.45561093090810223):0.37156673899045645):0.5356573388243997,27:0.5002436081418524):0.4382836447406966,(((66:0.009238168486191745,(83:1.097087549661701,59:0.035452835708333685):0.18557623088284902):0.06536010382226087,((87:0.9671241326044003,((16:0.5761335171481893,(50:0.08274841507119379,29:0.20452154014734525):0.008960605574878988):0.15466580443456246,6:0.031038985470871472):0.08873693148735529):0.4895952840353077,8:0.021237372146760913):0.029158925355997622):0.6265124178594199,20:0.0562256135785657):0.14155058845923874):0.34246719151547333,52:0.15238472576944284):0.4477326621284359,(81:0.21286177861903255,(38:0.42739746705523585,((32:0.2563939476177466,(77:0.09839158118983615,((46:0.3063079776955284,44:1.0764788865361645):0.0967881563934081,(64:0.015365891691861533,(37:0.3118178978180244,(18:1.0481894711877224,28:0.5696833775511019):0.10531476820290253):0.6750406605860493):0.2389474115301009):0.4460808727445591):0.16338377374941127):0.24658159668397372,(60:0.1941521696748576,(57:0.4396173353971071,45:0.14764439636767435):0.00292563495445064):0.5835606891253495):0.03559882517394586):0.007593162074185589):0.026737592907228436):0.5043586033493836):0.08892358714303672,(94:0.32167568499503174,13:0.24981928101855377):0.3005243433313902):1.9033152658984709,10:0.04014718449395094):0.7733848279809488):0.4247027874163072,26:0.028719487023320944):0.14531692560198994):0.46251150562929144,70:0.8772022746957981):0.20189500470051325):0.6201916487517867,(99:0.3985631513944936,(91:0.1098357173489104,(89:0.438021223588029,((((100:1.2155190484767608,((67:0.9245086601024375,61:0.9062831075622751):0.8799529224721387,(((19:0.6263320501489607,73:0.32725924441884124):0.545566863874281,((41:0.14479649843775055,76:0.35888793706597066):0.4751385103081651,(4:0.43409303833894364,97:0.03209835708235764):0.06955418065813235):0.20042318205952347):0.029130580497363834,75:0.4079664634631449):0.8342155556700948):0.22195875614900462):2.8313129468083824,((93:0.0032588777353970144,71:0.3314745510222519):0.18972134813675545,2:0.9790094658343982):0.2531977560327814):0.3212678043966397,((30:0.42125777935013353,(((92:0.23264294534966057,(((53:0.046927847524238536,1:0.12892356304959043):0.9182746307439067,((88:0.45401944305455244,12:0.07726464080803563):0.03531296007223317,24:0.3408916676707102):0.07312560201814489):0.8064335697234057,25:0.06574117610264962):0.3946280484435778):1.005385484039453,((80:0.10085302357579451,69:0.5381814163465224):0.18676540711800804,65:1.0653333275501993):0.22994806447027383):1.2907938991943295,(85:0.2786810201584169,43:0.3631279015020863):0.5591846264485989):0.370092595392598):0.01894215495227325,((((((82:0.07699009030211812,(11:0.019967389878505104,72:1.1986378804373903):0.7206619673334742):0.030554149050260992,96:0.06375870671885409):0.35439458784869826,(((86:0.5215929147365141,17:0.020839534194150744):0.22913355675960778,15:0.4899882008930847):0.038287250687130125,5:1.1227346249945698):0.04397277157006885):1.0508341848072957,48:0.2514002382673981):0.02272355663925163,(((((39:0.9008867900722497,90:0.5652518469346437):0.22912950939533494,((42:0.5288472434347464,((36:0.44411185594081104,55:0.6900228670372428):0.5030156268758961,23:0.9444877196027027):0.3582066070765175):0.06253855573049583,68:0.046825650993242096):0.056714246164616844):0.11292453687408499,(((62:0.1136023999488942,63:0.5282183994014511):0.08584662186109537,40:0.15052685144149636):0.09561535958851475,((35:1.0688079540107962,3:0.9323032634016541):0.02676955012696025,(51:0.9729500524306913,14:0.04486630061190611):0.1850717537389439):0.3980136894774997):0.4815062511338304):1.5066517480878021,31:0.24181583776313254):0.17640506748799334,(9:0.5292782083049792,21:0.06130431208511533):0.44107173512354025):0.11687497159870386):0.3801421398812672,33:0.0730333043861151):0.6999711958216572):0.5822070189933886):0.6830709433466415,(79:0.7470133542750927,54:1.495896451713568):0.0044491490776241704):0.2935533886021364):0.7140255334236354):0.9678607659741312):0.0675470728377956);",
//            /* 8 */  "(((((96:0.19593269703329863,(100:0.7544992453177741,24:0.032608961543205695):0.16801348050329512):0.217143210622309,(((80:0.001517617325462517,(((60:0.02600794849135557,((97:0.08862478248090877,(19:0.22693797531889626,14:1.1295478496913884):0.7510520098391376):0.5247595039550568,(((17:0.08377600437579646,((45:0.10555732175758381,((40:0.4671746347506913,58:1.1710806360458348):0.237996954915197,15:1.6658253812263526):0.9946316920992826):0.4687370662276278,(((((34:1.7724867001334168,49:1.4163660136774823):0.09799993689364817,((10:0.3009743358225263,(35:0.04773170379064973,(5:0.17885173014166966,16:0.7106681518198553):0.37256723224254795):0.053666752632955284):0.04099274904851291,(69:1.4439333283616942,((6:0.417897809280078,(71:0.02952484466050409,(54:0.10021888774073617,(59:0.43534404270662674,68:0.5665743072448134):0.06808936191880166):0.2418981434934132):0.14471019000705265):0.1285601271762129,((67:0.3519008103609327,36:0.3524636433882793):0.3665967423770632,23:0.20115734537071395):0.31480634945479935):0.209948045622246):0.5080298795131808):0.15986622454723554):0.6632809979575685,78:0.12343860143021823):0.2738378845452156,((27:0.0430531585652596,21:0.5598360365943802):0.35362454473229255,50:0.05979211997313172):0.24690107828472208):0.015925882324924423,(25:1.1164298102831889,29:0.93183326688396):0.618319524515369):0.06703753924466804):0.13366252535421363):0.5339325939095994,((((91:0.2772414826977889,52:0.5054410298692358):0.012245961880808842,(57:0.6908957477779047,46:0.1242660536287643):0.45678045518403665):0.09689181618779052,((43:0.05875674734537528,90:0.2470157776549904):0.6470927718147148,((((28:0.7062726345696688,41:0.12619466507616672):0.5279339445274762,((1:0.16197752379023367,12:0.0936878341766092):1.0118107676816173,4:0.06371547275534928):0.06331178642963398):0.1449167827806901,63:0.1231900161960271):0.23828417706345473,20:0.0378601589214822):0.045611678615212314):0.8587293083448599):0.5420674642468035,(64:0.12296240685548643,((26:1.3055554366821998,((22:0.7135388026823565,53:0.6553794907233433):0.5635353706486713,3:0.319169046517507):0.47913944032663514):0.6382409074552566,9:0.7129345589733305):0.28744783392818807):0.17306963200514502):0.7405187895273713):0.33399387891050125,8:1.8023038546554835):0.5924678803112702):0.5445164601831172):0.007238022980884828,32:1.3128883830865705):0.13340690953951384,(((98:0.523994867853685,((((81:0.13808210286308853,70:0.24449839687404396):0.09407430160994723,(39:0.02534540096873017,((89:0.20045498756714775,(83:0.08503372751564287,73:0.050892711823579906):0.14527900345950628):0.1007758060653059,87:0.6794603298260105):0.1486142740763876):0.9339744334947284):0.33286974301166516,(94:0.5872025248543964,(((99:1.130031886234744,65:0.0334368758185879):0.05493559297608952,33:1.4714983028992503):0.05406628655084811,11:0.4915430641248135):0.367449567894214):0.3933625310757476):0.7105188300285632,(66:0.07264775680667279,56:0.37051904950175896):0.08188191033372494):1.5059372514762073):0.2558116145927736,37:0.006350741961370865):0.09672071241581914,(77:0.19759715629071817,(79:0.12149768625576307,31:1.03643365865962):0.39015421362521874):0.3143643683653279):0.42068578176042415):0.09655626349881441):0.02460400823353126,86:0.36084437426554317):0.25306700347709654,72:0.82290814071923):0.6134198703925042):0.22465805069410472,82:0.09158964290059579):0.22023484997064813,(((30:0.017154279422356034,((((75:0.49143357922184805,(76:0.28537001523743943,(47:0.13980148118310265,88:0.27558089705798616):0.6356770970597386):0.7139308280815952):0.15999410869172204,((62:1.120571203782646,38:0.44549750204481153):0.15854167840146105,84:0.977150519578527):0.6993202553505358):0.5667700294334432,48:0.30455344416372787):2.3977573269791823,((((74:0.1257844219105282,(((13:0.26846285962892447,93:0.34994306207577575):0.0668117380329103,7:1.6134247531587778):0.06431121181186761,85:1.0596553228311065):0.36041129787105053):0.08853526300450731,18:0.7823801049245711):0.11583807527929313,42:0.1460761985367034):0.3036681077746568,2:0.6345246941343614):0.14627030205926506):0.0666910115268129):0.006447006748649287,95:0.33013306497314243):0.13748684668120426,44:0.9453930357994347):1.498888465149748):0.19473209209458542,(92:0.797016413269569,((51:0.7299446646702386,55:0.8426349822491348):0.2552910980781258,61:0.316238194437501):0.3443810869076849):1.3395752515186246);",
//            /* 9 */  "((72:0.26778410529745145,16:0.9099901060236721):0.572556639307777,((((95:0.16476383779159676,88:0.043010502026175956):0.027373721572081067,(49:0.5337319419158533,((65:0.7437555987541389,(((25:0.06520160789642127,2:0.09276120923381193):0.49176437927810923,((((66:0.1640592271541621,(30:0.6310067611534267,63:0.053353057470609144):0.26775377484213436):0.6911153103988585,23:0.25886630953279965):0.23672061440817593,85:0.6811271712242579):0.88047109253489,21:0.035509833182325146):0.1334075747734209):0.045841843378031655,(91:0.7810796396220372,((((35:0.3696497234316247,((99:0.2841562566009217,((((90:0.2091892549323302,93:0.39353395516714484):0.4232206217777461,9:0.11251447994286279):0.1800308940966806,((18:0.15316383890120688,69:0.14033738659421358):0.1844231030866611,(81:0.43578009647324467,73:0.3103644319971215):0.6485283802451252):1.069977767390587):0.7129782655656447,64:0.18449503942191292):0.4062398841393562):0.29321111632547225,48:0.0575819644497404):1.2387494936054075):0.3393774536470353,55:0.11323963232805934):0.07387562694960081,10:0.036610352098754184):0.03395775996906458,((((52:0.38144323157783244,80:0.7987611381333588):0.509346114812463,((((8:0.22586585677875637,53:0.3373759652024715):0.46455286573529553,((59:0.18862934983024987,1:0.4528645015033979):0.911575221871987,(56:0.5374669476738119,(17:0.3727190362752084,28:0.12812985951756062):0.17146441607895202):0.07484172326864602):0.526560218829603):0.18246989560362858,(44:1.763903077040106,(7:0.19339905954378656,20:0.3260783550961785):1.5124538821696814):0.08241196559786013):0.14021012116354425,58:0.5590182652099926):1.4398893973627445):0.34905785033321957,((74:0.09039135729219083,78:0.20995854952773785):0.5205312190594968,(24:0.3668850104763135,4:0.5145175498867807):0.7434984872360539):1.613008904745389):0.6102349600725585,(62:0.6044803001138535,((((((96:0.5072937801557419,31:0.5794214212731033):0.2712081265480051,68:0.8506805096120594):0.5974905834493063,(57:0.10700636565501898,51:0.5384611852745366):1.395093109724256):0.5019533182615326,(((((84:0.317204073514155,(((34:0.6175993994922453,32:0.7896059037064286):0.17947059708583135,(12:0.09115432415292446,60:0.47504664684058123):0.32326923070610825):0.42274952341622196,22:0.5197693871547759):0.08704082139413138):0.0088231898539044,(76:0.3938693591201038,94:0.37462781851914506):0.7773782014408085):0.020621354256771696,(61:1.3016030558283005,47:1.5796513156674639):0.04251393298459716):0.24836295995452384,89:0.2939297831039678):0.18415456033238842,(83:0.11559894501584678,82:0.5752882967091613):0.035998368754705634):0.36381139913404903):0.5307313823214246,36:0.43493768953607015):0.5686164880390647,(((((67:0.5287168805848876,(29:0.5197748572102281,(19:1.0258515788924416,41:0.10253767736666797):0.28420174331847314):0.04755892859287769):0.3028227898993232,13:0.9761816696651464):0.14466768486781945,43:0.08487172864726):0.12262855875809464,50:0.28590290495471304):0.1514980456284012,((92:0.04957813532763389,(37:0.008643063730506828,39:0.06388175496913862):0.216977791613008):0.07449478180585478,98:1.2060833459990103):0.6268705601455311):0.5501754535605508):0.014143451145258279):1.0077776094334494):0.6480814593866651):1.8610917060340562):0.23694373939904612):0.012096770663862166):1.1499279830077036,15:0.0989446231713007):0.09782579398140889):0.31599860360596743):0.8609861644022807,71:0.044818540117461225):0.05016986176329574,(((((((27:0.9192667546083975,(26:0.0900591654577454,46:0.19904964324240204):0.17486983722355198):0.03233681665933297,38:0.20545605518965893):0.3487562172361347,33:1.829218161398523):0.04487328698402493,(97:0.30044578766684804,87:0.050700523096797845):0.06054292682925411):0.2851509325367987,((((11:0.20798266945753596,(77:0.45443260432308286,5:0.6776501421086194):1.2052960699230568):0.7021089990669229,(100:0.036641296543967306,79:0.48749059631552927):2.0951639964127784):0.874795550652919,3:0.6104497577068013):0.5962534685036696,54:0.4761069361614636):0.5090942166401735):0.7364219385434367,6:0.2317442838428354):3.588516992657726,(((14:0.33085839693868735,42:0.1359545497645671):0.16976090503858643,(((75:0.3232483960272008,86:0.03351890552827985):0.12620180003638115,40:0.6959311697498141):0.4808677568002908,70:0.03828468578590982):0.29147282423511633):0.02294399302483452,45:1.9835733808193368E-4):1.400094368080568):0.568259308159158):0.19020027103967863);",
//            /* 10*/  "(((74:0.32981227091371323,(61:0.04503308082760071,((((97:0.4415529009725141,(91:0.9798103608717241,27:0.007745587863913883):0.07942990106969461):1.2769159284618148,28:0.018357012680464635):0.07421385884018328,67:0.010988614673328456):0.20269918955757404,((54:0.33424803276537496,(53:0.6833319504672986,(((((30:0.7149204751491034,19:0.1328998609477705):0.5286161357202213,96:1.192845854661313):1.5426411326305987,((26:0.46558316702624136,(71:1.2838507070550418,(4:1.8471165386594222,(((85:0.23372337345568406,(81:0.13666514053050038,5:0.6153505091584774):0.05903098010793073):0.855841826365686,83:0.31722231652923405):0.005665445735901997,(((((33:0.9875367622781843,99:0.5035332737307725):0.025989565009862936,84:0.8338312588209604):1.29743494765373,2:0.22374875572575625):0.9465989713141041,((((88:0.27236320760881805,11:0.11896248698604683):0.39506834236145805,(94:0.6073728052553622,42:0.31804837489583715):0.6207294375416619):0.04637568555454141,92:0.30800096620746364):0.08408374712585642,14:0.055289537060422056):0.4662809740077378):1.0548882541426607,46:0.3549566157230615):0.2989958819454621):0.10442819181779406):0.11673214955737166):0.27408086791743536):0.2875395174855493,(29:0.25219919926970213,(24:0.31075253085249077,25:0.6357825103706789):0.9347007949517891):0.22222181551988918):0.04092256999724153):1.2584507979295116,(70:0.019434479033757945,((44:0.0062545030920668765,((95:1.047397366664728,(52:0.34006529784617356,((((((47:1.0400425061198164,34:0.12282787949619411):0.3003582696154856,65:1.255615266857924):0.28061019974208445,(15:0.40991306473145883,41:2.5325759764251714):0.29883074949672706):0.08869241200262046,(49:0.35348462562931315,((48:0.5952630579869904,(72:0.4514206570177155,79:1.4283242213442708):0.19871948608119716):0.27548144452960877,12:0.892984722803791):0.44455045551456385):0.296860413392102):0.1641705162239746,((8:0.46686476673516797,((57:0.24179054456661264,10:0.010065770371418425):1.3076352901359078,(62:0.3049877998249979,((((((39:0.6593336279735633,80:0.16215091538111204):0.13671373823858168,((1:0.5374889016511959,98:0.3706747368135581):0.31640704372638795,64:0.5378967882260994):0.404209319450642):0.18200789367090042,90:0.9945839246994226):0.18799319544925197,45:0.35993176598759224):0.4805001706579122,(87:0.2299316509284952,20:1.8154332925124594):0.06488606618201231):0.19510829753371395,86:0.06894472708323773):0.3373684025347101):0.0056959492998269745):0.25231491181126486):0.12181032162614702,82:0.32401681418757766):0.44585600292315597):0.8964551643568837,(9:0.8182299720384876,(((16:0.7530769694312919,58:0.25659887324514496):0.8858547153564256,((13:0.7720021319127448,21:0.1809826867886699):0.16602736502837678,32:0.4048244114238031):0.7617183505627576):0.14137643667124378,(((((75:0.6503368892661168,36:0.6681993990037137):0.07734336586798518,(23:0.947923270388612,38:0.9574489826546935):0.06515522359532677):0.36911814314844116,77:0.2856180040191263):0.9187372332781449,((18:0.34868630512202475,50:0.26007250244097385):1.1632665967474327,3:1.4886189059882884):0.7676441391769049):0.016672870838133402,66:0.1689604150664259):0.14061842083130127):1.794793960408692):0.01642820701990111):0.8368305846625717):0.16123780567983914):0.035022606357221875,((89:2.3368911898605447,((76:0.649760192186692,63:0.01916583588484677):0.18707781548642233,(78:0.12717634239246678,31:0.043787753691652664):0.20409136472281997):0.17057872169254384):1.0136507337418088,(55:0.19514947568087493,17:0.13009194136081348):0.41479966211421093):0.24402253889293046):0.133397038420215):0.11414081547087918,(59:1.7109392994601142,((7:2.82622056805058,((68:0.11521280961551272,(37:0.3961142869085892,((22:0.09245304595081436,(93:0.10854203540491394,51:0.157205892369487):0.12374149365656104):0.16768907611151262,40:0.450066379049208):0.9217230804189102):0.05866460961760889):1.4357542152895868,69:0.12691297946788405):0.1762020854379247):0.7135474208651922,35:0.22156120885860586):0.3109217127373065):1.5386593468067389):0.9799666335924311):0.07876809613355107):0.6217290638352955,43:0.008058622432427498):2.1261353999875325):0.05438125770080404):0.398909900962483,100:0.31809473558773327):1.1788777106096084):0.24757224447711934):0.3981175105591799):0.7398874872692467,(6:0.03241667394145864,60:0.1638714537340764):0.628752104994776):0.1658371918015913,(73:0.1322536908343963,56:0.29379799240773075):0.9210859556235569);"
//    };
}
//        double sum = 0;
//        for (int n = 2; n <= 100; n++) {
//             sum += 1/(double) n;
//        }
//        System.out.println("sum = " + sum);

//    private static void getDate(int index, GetDateFromTree getDateFromTree, String treeString) { // single tree or import trees
//        Taxa taxa = new Taxa();
//        for (int n = 1; n <= 100; n++) {
//            Taxon t = new Taxon(Integer.toString(n));
//            taxa.addTaxon(t);
//        }
//        double[] tips = new double[taxa.getTaxonCount() + 1];
//        double rootHeight;
//        try {
//            Tree tree = getDateFromTree.importTree(taxa);
//            System.out.println(tree);
//
//            for (int i = 0; i < tree.getTaxonCount(); i++) {
//                FlexibleNode node = (FlexibleNode) tree.getExternalNode(i);
////                System.out.println(node.getTaxon() + " has " + node.getHeight());
//                tips[Integer.parseInt(node.getTaxon().getId())] = node.getHeight();
//            }
//            rootHeight = ((FlexibleNode) tree.getRoot()).getHeight();
//            System.out.println("tree " + index + " root height = " + rootHeight);
//            System.out.println("\n");
//
//        } catch (ImportException e) {
//            System.err.println("Error Parsing Input Tree: " + e.getMessage());
//            return;
//        } catch (IOException e) {
//            System.err.println("Error Parsing Input Tree: " + e.getMessage());
//            return;
//        }
//
//        if (index < 0) {
//            printXML(tips);
//        } else {
//            DecimalFormat twoDForm = new DecimalFormat("#0.##");
//            try {
//                outputXML("", index, tips, Double.valueOf(twoDForm.format(rootHeight + 1.0)), treeString);
//            } catch (IOException e) {
//                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
//        }
//        System.out.println("\n");
//    }
