/*
 * DebugUtils.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.mcmc;

/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @author Guy Baele
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

import dr.app.tools.TreeAnnotator;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.BranchRates;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.tree.TreeDebugLogger;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.loggers.Logger;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.math.MathUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class DebugUtils {

    private static final boolean DEBUG = false;

    /**
     * Writes out the current state in a human readable format to help debugging.
     * If it fails, then returns false but does not stop.
     * @param file the file
     * @param loggers loggers being used by MCMC of which we can use information to write trait-annotated trees to file
     * @param state the current state number
     * @param operatorSchedule
     * @return success
     */
    public static boolean writeStateToFile(File file, Logger[] loggers, long state, double lnL, OperatorSchedule operatorSchedule) {
        OutputStream fileOut = null;
        try {
            fileOut = new FileOutputStream(file);
            PrintStream out = new PrintStream(fileOut);

            int[] rngState = MathUtils.getRandomState();
            out.print("rng");
            for (int i = 0; i < rngState.length; i++) {
                out.print("\t");
                out.print(rngState[i]);
            }
            out.println();

            out.print("state\t");
            out.println(state);

            out.print("lnL\t");
            out.println(lnL);

            for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {
                out.print(parameter.getParameterName());
                out.print("\t");
                out.print(parameter.getDimension());
                for (int dim = 0; dim < parameter.getDimension(); dim++) {
                    out.print("\t");
                    out.print(parameter.getParameterValue(dim));
                }
                out.println();
            }

            for (int i = 0; i < operatorSchedule.getOperatorCount(); i++) {
                MCMCOperator operator = operatorSchedule.getOperator(i);
                out.print(operator.getOperatorName());
                out.print("\t");
                out.print(operator.getAcceptCount());
                out.print("\t");
                out.print(operator.getRejectCount());
                if (operator instanceof CoercableMCMCOperator) {
                    out.print("\t");
                    out.print(((CoercableMCMCOperator)operator).getCoercableParameter());
                }
                out.println();
            }

            for (Model model : Model.CONNECTED_MODEL_SET) {
                if (model instanceof TreeModel) {
                    out.print(model.getModelName());
                    out.print("\t");

                    //TODO This should be commented out at some point
                    out.println(((TreeModel) model).getNewick());

                    //make sure each TreeDebugLogger writes to a different file
                    for (Logger treeOutput : loggers) {
                        if (treeOutput instanceof TreeDebugLogger) {
                            ((TreeDebugLogger) treeOutput).writeToFile(state);
                        }
                    }
                }

            }

            out.close();
            fileOut.close();
        } catch (IOException ioe) {
            System.err.println("Unable to write file: " + ioe.getMessage());
            return false;
        }

        if (DEBUG) {
            for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
                System.err.println(likelihood.getId() + ": " + likelihood.getLogLikelihood());
            }
        }

        return true;
    }

    /**
     * Attempts to read the current state from a state dump file. This should be a state
     * dump created using the same XML file (some rudimentary checking of this is done).
     * If it fails then it will throw a RuntimeException. If successful it will return the
     * current state number.
     * @param file the file
     * @param loggers collection of loggers associated with the XML and debugging
     * @return the state number
     */
    public static long readStateFromFile(File file, Logger[] loggers, OperatorSchedule operatorSchedule, double[] lnL) {
        long state = -1;

        TreeParameterModel backup = null;

        try {
            FileReader fileIn = new FileReader(file);
            BufferedReader in = new BufferedReader(fileIn);

            int[] rngState = null;

            String line = in.readLine();
            String[] fields = line.split("\t");
            if (fields[0].equals("rng")) {
                // if there is a random number generator state present then load it...
                try {
                    rngState = new int[fields.length - 1];
                    for (int i = 0; i < rngState.length; i++) {
                        rngState[i] = Integer.parseInt(fields[i + 1]);
                    }

                } catch (NumberFormatException nfe) {
                    throw new RuntimeException("Unable to read state number from state file");
                }

                line = in.readLine();
                fields = line.split("\t");
            }

            try {
                if (!fields[0].equals("state")) {
                    throw new RuntimeException("Unable to read state number from state file");
                }
                state = Long.parseLong(fields[1]);
            } catch (NumberFormatException nfe) {
                throw new RuntimeException("Unable to read state number from state file");
            }

            line = in.readLine();
            fields = line.split("\t");
            try {
                if (!fields[0].equals("lnL")) {
                    throw new RuntimeException("Unable to read lnL from state file");
                }
                if (lnL != null) {
                    lnL[0] = Double.parseDouble(fields[1]);
                }
            } catch (NumberFormatException nfe) {
                throw new RuntimeException("Unable to read lnL from state file");
            }

            for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {

                line = in.readLine();
                fields = line.split("\t");
//                if (!fields[0].equals(parameter.getParameterName())) {
//                    System.err.println("Unable to match state parameter: " + fields[0] + ", expecting " + parameter.getParameterName());
//                }
                int dimension = Integer.parseInt(fields[1]);

                if (dimension != parameter.getDimension()) {
                    System.err.println("Unable to match state parameter dimension: " + dimension + ", expecting " + parameter.getDimension() + " for parameter: " + parameter.getParameterName());
                    System.err.print("Read from file: ");
                    for (int i = 0; i < fields.length; i++) {
                        System.err.print(fields[i] + "\t");
                    }
                    System.err.println();
                }

                if (fields[0].equals("branchRates.categories.rootNodeNumber")) {
                    // System.out.println("eek");
                    double value = Double.parseDouble(fields[2]);
                    parameter.setParameterValue(0, value);
                    if (DEBUG) {
                        System.out.println("restoring " + fields[0] + " with value " + value);
                    }
                } else {
                    if (DEBUG) {
                        System.out.print("restoring " + fields[0] + " with values ");
                    }
                    for (int dim = 0; dim < parameter.getDimension(); dim++) {
                        parameter.setParameterValue(dim, Double.parseDouble(fields[dim + 2]));
                        if (DEBUG) {
                            System.out.print(Double.parseDouble(fields[dim + 2]) + " ");
                        }
                    }
                    if (DEBUG) {
                        System.out.println();
                    }
                }

            }

            for (int i = 0; i < operatorSchedule.getOperatorCount(); i++) {
                MCMCOperator operator = operatorSchedule.getOperator(i);
                line = in.readLine();
                fields = line.split("\t");
                if (!fields[0].equals(operator.getOperatorName())) {
                    throw new RuntimeException("Unable to match operator: " + fields[0]);
                }
                if (fields.length < 3) {
                    throw new RuntimeException("Operator missing values: " + fields[0]);
                }
                operator.setAcceptCount(Integer.parseInt(fields[1]));
                operator.setRejectCount(Integer.parseInt(fields[2]));
                if (operator instanceof CoercableMCMCOperator) {
                    if (fields.length != 4) {
                        throw new RuntimeException("Coercable operator missing parameter: " + fields[0]);
                    }
                    ((CoercableMCMCOperator)operator).setCoercableParameter(Double.parseDouble(fields[3]));
                }
            }

            // load the tree models last as we get the node heights from the tree (not the parameters which
            // which may not be associated with the right node
            Set<String> expectedTreeModelNames = new HashSet<String>();
            for (Model model : Model.CONNECTED_MODEL_SET) {
                if (model instanceof TreeModel) {
                    if (DEBUG) {
                        System.out.println("model " + model.getModelName());
                    }
                    expectedTreeModelNames.add(model.getModelName());
                    if (DEBUG) {
                        for (String s : expectedTreeModelNames) {
                            System.out.println(s);
                        }
                    }
                }
            }

            int counter = 0;
            //TreeDebugLoggers are not created when loading a dump file so no way to get the file names
            //Get the file names using the File that is passed as an argument; may be unsafe
            ArrayList<String> fileNames = new ArrayList<String>();
            for (Logger logger : loggers) {
                if (logger instanceof TreeLogger) {
                    fileNames.add(file.getName() + ".tree." + counter);
                    counter++;
                }
            }

            counter = 0;

            // Read in all (possibly more than one) tree
            while((line = in.readLine()) != null) {
                fields = line.split("\t");
                boolean treeFound = false;

                for (Model model : Model.CONNECTED_MODEL_SET) {

                    if (DEBUG) {
                        System.out.println(model);
                    }

                    //TODO final line of dump file can be removed; keep it in for now but don't use it
                    if (model instanceof TreeModel && fields[0].equals(model.getModelName())) {
                        if (DEBUG) {
                            System.out.println("TREEMODEL");
                        }
                        treeFound = true;

                        TreeModel checkTree = ((TreeModel) model);

                        Set<String> attributeNames = new HashSet<String>();

                        if (DEBUG) {
                            for (int i = 0; i < checkTree.getNodeCount(); i++) {
                                NodeRef node = checkTree.getNode(i);
                                System.out.print(node + ": ");
                                Iterator iter = checkTree.getNodeAttributeNames(node);
                                if (iter != null) {
                                    while (iter.hasNext()) {
                                        String name = (String) iter.next();
                                        System.out.print(name + " ");
                                        attributeNames.add(name);
                                    }
                                }
                                System.out.println();
                            }
                        }

                        if (backup != null && DEBUG) {
                            System.out.println("Checking rate category assignments in TreeParameterModel before adopting tree structure:");
                            for (int i = 0; i < checkTree.getNodeCount(); i++) {
                                if (checkTree.getParent(checkTree.getNode(i)) == null) {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- null");
                                } else {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- " + checkTree.getParent(checkTree.getNode(i)).getNumber() + ": " + backup.getNodeValue(checkTree, checkTree.getNode(i)));
                                }
                            }
                        }

                        System.out.println("Attempting to read tree file: " + fileNames.get(counter));

                        //boolean needs to be set to false to read node attributes
                        TreeImporter importer = new NexusImporter(new FileReader(fileNames.get(counter)), false);
                        Tree tree = importer.importNextTree();

                        counter++;

                        if (DEBUG) {
                            System.out.println(tree);
                        }

                        if (importer.hasTree()) {
                            throw new RuntimeException("DebugUtils: only one tree can be loaded from a saved file.");
                        }


                        if (DEBUG) {
                            attributeNames = new HashSet<String>();

                            for (int i = 0; i < tree.getNodeCount(); i++) {
                                NodeRef node = tree.getNode(i);
                                System.out.print(node + ": ");
                                Iterator iter = tree.getNodeAttributeNames(node);
                                if (iter != null) {
                                    while (iter.hasNext()) {
                                        String name = (String) iter.next();
                                        System.out.print(name + " ");
                                        attributeNames.add(name);
                                    }
                                }
                                System.out.println();
                            }

                            for (int i = 0; i < tree.getNodeCount(); i++) {
                                System.out.print(tree.getNode(i).getNumber() + " (" + tree.getNodeTaxon(tree.getNode(i)) + ") " + " -- ");
                                if (!tree.isRoot(tree.getNode(i))) {
                                    System.out.print(tree.getParent(tree.getNode(i)).getNumber() + ": ");
                                }
                                for (String attributeName : attributeNames) {
                                    Object value = tree.getNodeAttribute(tree.getNode(i), attributeName);
                                    System.out.print(value + " ");
                                }
                                System.out.println();
                            }

                        }

                        //adopt the loaded tree structure; this does not yet copy the traits on the branches
                        ((TreeModel) model).beginTreeEdit();
                        ((TreeModel) model).adoptTreeStructure(tree);
                        ((TreeModel) model).endTreeEdit();

                        expectedTreeModelNames.remove(model.getModelName());

                        //TODO Still set the rate categories to each Node
                        //TODO But first check them in the TreeModel (i.e. not in the Tree)

                        if (backup != null && DEBUG) {
                            System.out.println("Checking rate category assignments in TreeParameterModel after adopting tree structure:");
                            for (int i = 0; i < checkTree.getNodeCount(); i++) {
                                if (checkTree.getParent(checkTree.getNode(i)) == null) {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- null");
                                } else {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- " + checkTree.getParent(checkTree.getNode(i)).getNumber() + ": " + backup.getNodeValue(checkTree, checkTree.getNode(i)));
                                }
                            }
                        }

                        //If we have found a TreeParameterModel before, we'll find it again here
                        //TODO Set the rates according to those in the donor tree!
                        if (backup != null) {

                            if (DEBUG) {
                                System.out.println("Copying rate categories:");
                            }

                            for (int i = 0; i < checkTree.getExternalNodeCount(); i++) {
                                //TODO Use method that takes a taxon and moves up to the root and set the rate categories along the way
                                //TODO Hence only iterate over external nodes!
                                //first find corresponding Taxon in donor tree
                                Taxon donorTaxon = getDonorTaxon(checkTree.getNodeTaxon(checkTree.getExternalNode(i)), tree);
                                //then set the attributes (for now: rate categories) from the tips up until the root
                                //TODO DiscretizedBranchRates.RATECATEGORY is currently hard-coded as an argument
                                //TODO Figure out other examples of TreeTraits for which this code can be used and get the traits names automatically
                                adoptAttributeValueFromDonorTree(checkTree, backup, tree, checkTree.getNodeTaxon(checkTree.getExternalNode(i)), donorTaxon, DiscretizedBranchRates.RATECATEGORY);
                            }

                        }

                        if (backup != null && DEBUG) {
                            System.out.println("Checking rate category assignments in TreeParameterModel after adopting rate categories:");
                            for (int i = 0; i < checkTree.getNodeCount(); i++) {
                                if (checkTree.getParent(checkTree.getNode(i)) == null) {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- null");
                                } else {
                                    System.out.println("  " + checkTree.getNode(i).getNumber() + " (" + checkTree.getNodeTaxon(checkTree.getNode(i)) + ") " + " -- " + checkTree.getParent(checkTree.getNode(i)).getNumber() + ": " + backup.getNodeValue(checkTree, checkTree.getNode(i)));
                                }
                            }
                        }

                    }

                    if (model instanceof TreeParameterModel) {
                        backup = ((TreeParameterModel) model);

                        //back up the TreeParameterModel for use in restoring the TreeModel
                        TreeModel temp = ((TreeParameterModel) model).getTreeModel();

                        if (DEBUG) {
                            System.out.println("TREEPARAMETERMODEL");
                            System.out.println("Node heights may be messed up from already setting the node heights previously");
                            System.out.println("node count: " + temp.getNodeCount());
                            System.out.println("root: " + temp.getRoot());

                            for (int i = 0; i < temp.getNodeCount(); i++) {
                                if (i == (temp.getNodeCount() - 1)) {
                                    System.out.println("  " + temp.getNode(i).getNumber() + " (" + temp.getNodeTaxon(temp.getNode(i)) + ") " + " -- null");
                                } else {
                                    System.out.println("  " + temp.getNode(i).getNumber() + " (" + temp.getNodeTaxon(temp.getNode(i)) + ") " + " -- " + temp.getParent(temp.getNode(i)).getNumber() + ": " + backup.getNodeValue(temp, temp.getNode(i)));
                                }
                            }
                        }

                    }

                }

                if (!treeFound) {
                    throw new RuntimeException("Unable to match state parameter: " + fields[0]);
                }
            }

            if (expectedTreeModelNames.size() > 0) {
                StringBuilder sb = new StringBuilder();
                for (String notFoundName : expectedTreeModelNames) {
                    sb.append("Expecting, but unable to match state parameter:" + notFoundName + "\n");
                }
                throw new RuntimeException(sb.toString());
            }

            if (DEBUG) {
                System.out.println("\nDouble checking:");
                for (Parameter parameter : Parameter.CONNECTED_PARAMETER_SET) {
                    if (parameter.getParameterName().equals("branchRates.categories.rootNodeNumber")) {
                        System.out.println(parameter.getParameterName() + ": " + parameter.getParameterValue(0));
                    }
                }
            }

            if (rngState != null) {
                MathUtils.setRandomState(rngState);
            }

            in.close();
            fileIn.close();
            for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
                likelihood.makeDirty();
            }
        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
        } catch (Importer.ImportException ie) {
            throw new RuntimeException("Unable to import tree: " + ie.getMessage());
        }

        return state;
    }

    /**
     *
     * @param recipient The recipient tree to which the attribute value must be copied
     * @param tpm TreeParameterModel corresponding to the recipient TreeModel
     * @param donor The donor tree from which to copy the attribute value
     */
    private static void adoptAttributeValueFromDonorTree(TreeModel recipient, TreeParameterModel tpm, Tree donor, Taxon recipientTaxon, Taxon donorTaxon, String attributeName) {

        if (DEBUG) {
            System.out.println(recipientTaxon.toString());
            System.out.println(donorTaxon.toString());
        }

        NodeRef recipientNode = null;
        NodeRef donorNode = null;
        //TODO Don't like this but works for now
        for (int i = 0; i < recipient.getExternalNodeCount(); i++) {
            if (recipient.getNodeTaxon(recipient.getExternalNode(i)).equals(recipientTaxon)) {
                recipientNode = recipient.getExternalNode(i);
                break;
            }
        }
        for (int i = 0; i < donor.getExternalNodeCount(); i++) {
            if (donor.getNodeTaxon(donor.getExternalNode(i)).equals(donorTaxon)) {
                donorNode = donor.getExternalNode(i);
                break;
            }
        }

        if (DEBUG) {
            System.out.println(recipientNode.getNumber());
            System.out.println(donorNode.getNumber());
        }

        //TODO check if we need this to work for non-numerical attributes
        while (donor.getParent(donorNode) != null) {
            //System.out.println(tpm.getNodeValue(recipient, recipientNode));
            tpm.setNodeValue(recipient, recipientNode, (Integer)donor.getNodeAttribute(donorNode, attributeName));
            //System.out.println(recipientNode + ": " + (Integer)donor.getNodeAttribute(donorNode, attributeName));
            //System.out.println(tpm.getNodeValue(recipient, recipientNode));
            recipientNode = recipient.getParent(recipientNode);
            donorNode = donor.getParent(donorNode);
        }

    }

    private static Taxon getDonorTaxon(Taxon taxon, Tree tree) {

        //System.out.println("Searching for taxon: " + taxon);
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            if (taxon.getId().equals(tree.getNodeTaxon(tree.getNode(i)).getId())) {
                return tree.getNodeTaxon(tree.getNode(i));
            }
        }

        //if a taxon is not found, throw an exception for now
        //TODO Not needed when number of taxa is growing
        throw new RuntimeException("Taxon " + taxon.toString() + " not found!");

    }


}
