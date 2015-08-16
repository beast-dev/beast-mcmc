/*
 * GenerateRelaxedClockXMLByData.java
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

import dr.app.beauti.util.XMLWriter;

import java.io.*;

/**
 * Given data file format is
 * 1st line : anything which is ignored
 * 2nd line : one taxon name    separator (e.g. space or tab)    its sequence data
 * 3rd line : ...
 *
 * @author Walter Xie
 */
public class GenerateRelaxedClockXMLByData {

    static final String path = "/Users/local/EC/dxie004/Documents/BEAST1Release/RR_Datasets/";
    static final String startingTree = "(((((Ssc:65,Bta:65):16,((Cfa:46,Fca:46):28,Eca:74):7):11," +
            "(((Rno:20,Mmu:20):65,Ocu:85):5,(((Hsa:5,Ptr:5):5,Ppy:10):13,Mml:23):67):2):81,Tvu1:173):137,Gga:310);";
    static final double rescaleHeight = 0.5;

    static public void main(String[] args) {
        int taxonNum = 14;

        for (int c = 1; c <= 100; c++) {

            String inputFileName = "geneRRNP" + c + ".phy";
            String[][] data = new String[taxonNum][2];

            try {
                System.out.println("Input data from " + path + inputFileName + "\n\n");
                FileReader fileReader = new FileReader(path + inputFileName);

                LineNumberReader lineNumberReader = new LineNumberReader(fileReader);

                String line = "";
                int lineNum = 0;
                while ((line = lineNumberReader.readLine()) != null) {
                    lineNum = lineNumberReader.getLineNumber();
                    System.out.println("Line:  " + lineNum + ": " + line);

                    if (lineNumberReader.getLineNumber() > 1) {  // not need 1st row

                        String[] d = line.split("\\s+");

                        if (d.length != 2) throw new IOException("d " + d.length);

                        data[lineNum - 2][0] = d[0]; // taxon name
                        data[lineNum - 2][1] = d[1]; // alignment
                    }
                }
                if (lineNum != 15) throw new IOException("not have 15 line " + lineNum);

                fileReader.close();

                String outputFileName = "geneRRNP" + c + ".xml";
                System.out.println("\n\nCreating xml : " + path + outputFileName);
                XMLWriter w = new XMLWriter(new BufferedWriter(new FileWriter(new File(path + outputFileName))));

                writeData(w, data);

                writeRestXML(w, outputFileName);

                w.close();
                System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++\n\n");

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
    }


    private static void writeRestXML(XMLWriter w, String outputFileName) throws IOException {

        w.writeText("\n\t<!-- a Birth-Death speciation process (Gernhard 2008).                       -->\n" +
                "\t<birthDeathModel id=\"birthDeath\" units=\"substitutions\">\n" +
                "\t\t<birthMinusDeathRate>\n" +
                "\t\t\t<parameter id=\"birthDiffRate\" value=\"1.0\" lower=\"0.0\" upper=\"1000000.0\"/>\n" +
                "\t\t</birthMinusDeathRate>\n" +
                "\t\t<relativeDeathRate>\n" +
                "\t\t\t<parameter id=\"relativeDeathRate\" value=\"0.5\" lower=\"0.0\" upper=\"1.0\"/>\n" +
                "\t\t</relativeDeathRate>\n" +
                "\t</birthDeathModel>\n");

//        w.writeText("\n\t<constantSize id=\"constant\" units=\"substitutions\">\n" +
//                "\t\t<populationSize>\n" +
//                "\t\t\t<parameter id=\"popSize\" value=\"0.077\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
//                "\t\t</populationSize>\n" +
//                "\t</constantSize>\n");

        w.flush();
        w.writeText("\n" +
                "\t<!-- Generate a random starting tree under the coalescent process      -->\n" +
                "\t<newick id=\"startingTree\" rescaleHeight=\"" + rescaleHeight + "\">\n");
        w.writeText(startingTree);
        w.writeText("\n" + "\t</newick>\n");

//        w.writeText("\n<!-- Construct a rough-and-ready UPGMA tree as an starting tree              -->\n" +
//                "\t<upgmaTree id=\"startingTree\">\n" +
//                "\t\t<distanceMatrix correction=\"JC\">\n" +
//                "\t\t\t<patterns>\n"+
//                "\t\t\t\t<alignment idref=\"alignment\"/>\n" +
//                "\t\t\t</patterns>\n" +
//                "\t\t</distanceMatrix>\n" +
//                "\t</upgmaTree>\n");

        w.flush();

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
                "\t</treeModel>\n");

        w.flush();
//        w.writeText("\n<!-- Generate a coalescent likelihood                                        -->\n" +
//                "\t<coalescentLikelihood id=\"coalescent\">\n" +
//                "\t\t<model>\n" +
//                "\t\t\t<constantSize idref=\"constant\"/>\n" +
//                "\t\t</model>\n" +
//                "\t\t<populationTree>\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t</populationTree>\n" +
//                "\t</coalescentLikelihood>\n");

        w.writeText("\n" +
                "\t<!-- Generate a speciation likelihood for Yule or Birth Death                -->\n" +
                "\t<speciationLikelihood id=\"speciation\">\n" +
                "\t\t<model>\n" +
                "\t\t\t<birthDeathModel idref=\"birthDeath\"/>\n" +
                "\t\t</model>\n" +
                "\t\t<speciesTree>\n" +
                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t</speciesTree>\n" +
                "\t</speciationLikelihood>\n");

        w.writeText("\n<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut (2006) PLoS Biology 4, e88 )-->\n" +
                "\t<discretizedBranchRates id=\"branchRates\">\n" +
                "\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t<distribution>\n" +
                "\t\t\t<logNormalDistributionModel meanInRealSpace=\"true\">\n" +
                "\t\t\t\t<mean>\n" +
                "\t\t\t\t\t<parameter id=\"ucld.mean\" value=\"1.0\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
                "\t\t\t\t</mean>\n" +
                "\t\t\t\t<stdev>\n" +
                "\t\t\t\t\t<parameter id=\"ucld.stdev\" value=\"0.3333333333333333\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
                "\t\t\t\t</stdev>\n" +
                "\t\t\t</logNormalDistributionModel>\n" +
                "\t\t</distribution>\n" +
                "\t\t<rateCategories>\n" +
                "\t\t\t<parameter id=\"branchRates.categories\" dimension=\"10\"/>\n" +
                "\t\t</rateCategories>\n" +
                "\t</discretizedBranchRates>\n");

//        w.writeText("\n\t<rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">\n" +
//                "\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t<discretizedBranchRates idref=\"branchRates\"/>\n" +
//                "\t</rateStatistic>\n" +
//                "\t<rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">\n" +
//                "\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t<discretizedBranchRates idref=\"branchRates\"/>\n" +
//                "\t</rateStatistic>\n" +
//                "\t<rateCovarianceStatistic id=\"covariance\" name=\"covariance\">\n" +
//                "\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t<discretizedBranchRates idref=\"branchRates\"/>\n" +
//                "\t</rateCovarianceStatistic>\n");


        w.writeText("\t<!-- The HKY substitution model (Hasegawa, Kishino & Yano, 1985)             -->\n" +
                "\t<HKYModel id=\"hky\">\n" +
                "\t\t<frequencies>\n" +
                "\t\t\t<frequencyModel dataType=\"nucleotide\">\n" +
                "\t\t\t\t<frequencies>\n" +
                "\t\t\t\t\t<parameter id=\"hky.frequencies\" value=\"0.25 0.25 0.25 0.25\"/>\n" +
                "\t\t\t\t</frequencies>\n" +
                "\t\t\t</frequencyModel>\n" +
                "\t\t</frequencies>\n" +
                "\t\t<kappa>\n" +
                "\t\t\t<parameter id=\"hky.kappa\" value=\"2.0\" lower=\"0.0\" upper=\"Infinity\"/>\n" +
                "\t\t</kappa>\n" +
                "\t</HKYModel>\n" +
                "\n" +
                "\t<!-- site model                                                              -->\n" +
                "\t<siteModel id=\"siteModel\">\n" +
                "\t\t<substitutionModel>\n" +
                "\t\t\t<HKYModel idref=\"hky\"/>\n" +
                "\t\t</substitutionModel>\n" +
                "\t</siteModel>\n");

        w.writeText("\t<!-- Likelihood for tree given sequence data                                 -->\n" +
                "\t<treeLikelihood id=\"treeLikelihood\" useAmbiguities=\"false\">\n" +
                "\t\t<patterns idref=\"patterns\"/>\n" +
                "\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t<siteModel idref=\"siteModel\"/>\n" +
                "\t\t<discretizedBranchRates idref=\"branchRates\"/> \n" +
                "\t</treeLikelihood>\n");

        w.flush();

        w.writeText("\n" + "\t<!-- Define operators                                                        -->\n" +
                "\t<operators id=\"operators\">\n");

        w.writeText("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n" +
                "\t\t\t<parameter idref=\"hky.kappa\"/>\n" +
                "\t\t</scaleOperator>\n" +
                "\t\t<deltaExchange delta=\"0.01\" weight=\"0.4\">\n" +
                "\t\t\t<parameter idref=\"hky.frequencies\"/>\n" +
                "\t\t</deltaExchange>\n" +
//                "\t\t<subtreeSlide size=\"0.0077\" gaussian=\"true\" weight=\"15\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t</subtreeSlide>\n" +
//                "\t\t<narrowExchange weight=\"15\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t</narrowExchange>\n" +
//                "\t\t<wideExchange weight=\"3\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t</wideExchange>\n" +
//                "\t\t<wilsonBalding weight=\"3\">\n" +
//                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
//                "\t\t</wilsonBalding>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"1\">\n" +
                "\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t</scaleOperator>\n" +
                "\t\t<uniformOperator weight=\"12\">\n" +
                "\t\t\t<parameter idref=\"treeModel.internalNodeHeights\"/>\n" +
                "\t\t</uniformOperator>\n" +
//                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n" +
//                "\t\t\t<parameter idref=\"popSize\"/>\n" +
//                "\t\t</scaleOperator>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n" +
                "\t\t\t<parameter idref=\"birthDiffRate\"/>\n" +
                "\t\t</scaleOperator>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">\n" +
                "\t\t\t<parameter idref=\"relativeDeathRate\"/>\n" +
                "\t\t</scaleOperator>\n" +

                "\t\t<upDownOperator scaleFactor=\"0.75\" weight=\"13\">\n" +
                "\t\t\t<up>\n" +
                "\t\t\t</up>\n" +
                "\t\t\t<down>\n" +
                "\t\t\t\t<parameter idref=\"treeModel.allInternalNodeHeights\"/>\n" +
                "\t\t\t</down>\n" +
                "\t\t</upDownOperator>\n" +
                "\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"3\">\n" +
                "\t\t\t<parameter idref=\"ucld.stdev\"/>\n" +
                "\t\t</scaleOperator> \n" +
                "\t\t<uniformIntegerOperator weight=\"13\">\n" +
                "\t\t\t<parameter idref=\"branchRates.categories\"/>\n" +
                "\t\t</uniformIntegerOperator>\n" +
                "\t\t<swapOperator size=\"1\" weight=\"13\" autoOptimize=\"false\">\n" +
                "\t\t\t<parameter idref=\"branchRates.categories\"/>\n" +
                "\t\t</swapOperator>\n");

        w.writeText("\t</operators>");

        w.flush();
        w.writeText("\n" +
                "\t<!-- Define MCMC                                                             -->\n" +
                "\t<mcmc id=\"mcmc\" chainLength=\"1000000\" autoOptimize=\"true\" " +
                "operatorAnalysis=\"" + outputFileName + ".ops\">\n" +
                "\t\t<posterior id=\"posterior\">\n" +
                "\t\t\t<prior id=\"prior\">\n" +
                "\t\t\t\t<logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">\n" +
                "\t\t\t\t\t<parameter idref=\"hky.kappa\"/>\n" +
                "\t\t\t\t</logNormalPrior>\n" +
                "\t\t\t\t<uniformPrior lower=\"0.0\" upper=\"1.0\">\n" +
                "\t\t\t\t\t<parameter idref=\"hky.frequencies\"/>\n" +
                "\t\t\t\t</uniformPrior>\n" +
                "\t\t\t\t<exponentialPrior mean=\"0.3333333333333333\" offset=\"0.0\">\n" +
                "\t\t\t\t\t<parameter idref=\"ucld.stdev\"/>\n" +
                "\t\t\t\t</exponentialPrior>  \n" +
//                "\t\t\t\t<oneOnXPrior>\n" +
//                "\t\t\t\t\t<parameter idref=\"popSize\"/>\n" +
//                "\t\t\t\t</oneOnXPrior> \n" +
//                "\t\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n");
                "\t\t\t<speciationLikelihood idref=\"speciation\"/>\n");

        w.writeText("\n" +
                "\t\t\t</prior>\n" +
                "\t\t\t<likelihood id=\"likelihood\">\n" +
                "\t\t\t\t<treeLikelihood idref=\"treeLikelihood\"/>\n" +
                "\t\t\t</likelihood>\n" +
                "\t\t</posterior>\n" +
                "\t\t<operators idref=\"operators\"/>\n");

        w.flush();

        w.writeText("\n" +
                "\t\t<!-- write log to screen                                                     -->\n" +
                "\t\t<log id=\"screenLog\" logEvery=\"100000\">\n" +
                "\t\t\t<column label=\"Posterior\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t</column>\n" +
//                "\t\t\t<column label=\"Prior\" dp=\"4\" width=\"12\">\n" +
//                "\t\t\t\t<prior idref=\"prior\"/>\n" +
//                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"Likelihood\" dp=\"4\" width=\"12\">\n" +
                "\t\t\t\t<likelihood idref=\"likelihood\"/>\n" +
                "\t\t\t</column>\n" +
                "\t\t\t<column label=\"rootHeight\" sf=\"6\" width=\"12\">\n" +
                "\t\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t</column>\n" +
//                "\t\t\t<parameter idref=\"ucld.mean\"/>\n" +
                "\t\t</log>\n"
        );


        w.writeText("\t\t<!-- write log to file                          -->\n" +
                "\t\t<log id=\"fileLog\" logEvery=\"100\" fileName=\"" + outputFileName + ".log\">\n" +
                "\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t\t<prior idref=\"prior\"/>\n" +
                "\t\t\t<treeLikelihood idref=\"treeLikelihood\"/>\n" +
//                "\t\t\t<coalescentLikelihood idref=\"coalescent\"/>\n" +
                "\t\t\t<speciationLikelihood idref=\"speciation\"/>\n" +
//                "\t\t\t<parameter idref=\"popSize\"/>\n" +
                "\t\t\t<parameter idref=\"birthDiffRate\"/>\n" +
                "\t\t\t<parameter idref=\"relativeDeathRate\"/>\n" +
                "\t\t\t<parameter idref=\"treeModel.rootHeight\"/>\n" +
                "\t\t\t<parameter idref=\"hky.kappa\"/>\n" +
                "\t\t\t<parameter idref=\"hky.frequencies\"/>\n" +
                "\t\t\t<parameter idref=\"ucld.mean\"/>\n" +
                "\t\t\t<parameter idref=\"ucld.stdev\"/>\n" +
                "\t\t\t<parameter idref=\"branchRates.categories\"/> \n");
        w.writeText("\t\t</log>\n");

        w.writeText("\t\t<logTree id=\"treeFileLog\" logEvery=\"10000\" nexusFormat=\"true\" " +
                "fileName=\"" + outputFileName + ".trees\" sortTranslationTable=\"true\">\n" +
                "\t\t\t<treeModel idref=\"treeModel\"/>\n" +
                "\t\t\t<discretizedBranchRates idref=\"branchRates\"/> \n" +
                "\t\t\t<posterior idref=\"posterior\"/>\n" +
                "\t\t</logTree>\n");

        w.writeText("\t</mcmc>\n" +
                "\t<report>\n" +
                "\t\t<property name=\"timer\">\n" +
                "\t\t\t<mcmc idref=\"mcmc\"/>\n" +
                "\t\t</property>\n" +
                "\t</report>\n" +
                "</beast>\n");

        w.flush();
    }

    private static void writeData(XMLWriter w, String[][] tips) throws IOException {
        w.writeText("<?xml version=\"1.0\" standalone=\"yes\"?>\n" + "\n" +
                "<!--       Generated by BEAUTi v1.6.2                                        -->\n" +
                "<!--       by Alexei J. Drummond and Andrew Rambaut                          -->\n" +
                "<!--       Department of Computer Science, University of Auckland and        -->\n" +
                "<!--       Institute of Evolutionary Biology, University of Edinburgh        -->\n" +
                "<!--       http://beast.bio.ed.ac.uk/                                        -->\n" +
                "<beast>\n" + "\n" +
                "\t<!-- The list of taxa to be analysed (can also include dates/ages).          -->\n" +
                "\t<!-- ntax=" + tips.length + "                                             -->\n");

        w.writeText("\t<taxa id=\"taxa\">\n");
        for (int n = 0; n < tips.length; n++) {
            w.writeText("\t\t<taxon id=\"" + tips[n][0] + "\"/>\n");
        }
        w.writeText("\t</taxa>\n");
        w.flush();

        w.writeText("\t<alignment id=\"alignment\" dataType=\"nucleotide\">\n");
        for (int n = 0; n < tips.length; n++) {
            w.writeText("\t\t<sequence>\n");
            w.writeText("\t\t<taxon idref=\"" + tips[n][0] + "\"/>\n");
            w.writeText(tips[n][1]);
            w.writeText("\t\t</sequence>\n");
            w.flush();
        }
        w.writeText("\t</alignment>\n");

        w.writeText("\n\t<patterns id=\"patterns\" from=\"1\">\n" +
                "\t\t<alignment idref=\"alignment\"/>\n" +
                "\t</patterns>");

    }

}
