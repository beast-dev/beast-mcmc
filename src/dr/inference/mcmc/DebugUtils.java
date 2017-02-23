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

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
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
import java.util.Set;

public class DebugUtils {

    private static final boolean DEBUG = false;

    /**
     * Writes out the current state in a human readable format to help debugging.
     * If it fails, then returns false but does not stop.
     * @param file the file
     * @param state the current state number
     * @param operatorSchedule
     * @return success
     */
    public static boolean writeStateToFile(File file, long state, double lnL, OperatorSchedule operatorSchedule) {
        OutputStream fileOut = null;
        try {
            fileOut = new FileOutputStream(file);
            PrintStream out = new PrintStream(fileOut);

            ArrayList<TreeParameterModel> traitModels = new ArrayList<TreeParameterModel>();

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
                out.print("parameter");
                out.print("\t");
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
                out.print("operator");
                out.print("\t");
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
                    out.print("tree");
                    out.print("\t");
                    out.println(model.getModelName());

                    //replace Newick format by printing general graph structure
                    //out.println(((TreeModel) model).getNewick());

                    out.println("#node height taxon");
                    int nodeCount = ((TreeModel) model).getNodeCount();
                    out.println(nodeCount);
                    for (int i = 0; i < nodeCount; i++) {
                        out.print(((TreeModel) model).getNode(i).getNumber());
                        out.print("\t");
                        out.print(((TreeModel) model).getNodeHeight(((TreeModel) model).getNode(i)));
                        if (((TreeModel) model).isExternal(((TreeModel) model).getNode(i))) {
                            out.print("\t");
                            out.print(((TreeModel) model).getNodeTaxon(((TreeModel) model).getNode(i)).getId());
                        }
                        out.println();
                    }

                    out.println("#edges");
                    out.println("#child-node parent-node");

                    out.println(nodeCount);
                    for (int i = 0; i < nodeCount; i++) {
                        if (((TreeModel) model).getParent(((TreeModel) model).getNode(i)) != null) {
                            out.print(((TreeModel) model).getNode(i).getNumber());
                            out.print("\t");
                            out.print(((TreeModel) model).getParent(((TreeModel) model).getNode(i)).getNumber());

                            for (TreeParameterModel tpm : traitModels) {
                                out.print("\t");
                                out.println(tpm.getNodeValue((TreeModel)model, ((TreeModel) model).getNode(i)));
                            }
                        }
                    }

                }

                if (model instanceof TreeParameterModel) {
                    traitModels.add((TreeParameterModel)model);
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
     * @return the state number
     */
    public static long readStateFromFile(File file, OperatorSchedule operatorSchedule, double[] lnL) {
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
                //if (!fields[0].equals(parameter.getParameterName())) {
                //  System.err.println("Unable to match state parameter: " + fields[0] + ", expecting " + parameter.getParameterName());
                //}
                int dimension = Integer.parseInt(fields[2]);

                if (dimension != parameter.getDimension()) {
                    System.err.println("Unable to match state parameter dimension: " + dimension + ", expecting " + parameter.getDimension() + " for parameter: " + parameter.getParameterName());
                    System.err.print("Read from file: ");
                    for (int i = 0; i < fields.length; i++) {
                        System.err.print(fields[i] + "\t");
                    }
                    System.err.println();
                }

                if (fields[1].equals("branchRates.categories.rootNodeNumber")) {
                    // System.out.println("eek");
                    double value = Double.parseDouble(fields[3]);
                    parameter.setParameterValue(0, value);
                    if (DEBUG) {
                        System.out.println("restoring " + fields[1] + " with value " + value);
                    }
                } else {
                    if (DEBUG) {
                        System.out.print("restoring " + fields[1] + " with values ");
                    }
                    for (int dim = 0; dim < parameter.getDimension(); dim++) {
                        parameter.setParameterValue(dim, Double.parseDouble(fields[dim + 3]));
                        if (DEBUG) {
                            System.out.print(Double.parseDouble(fields[dim + 3]) + " ");
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
                if (!fields[1].equals(operator.getOperatorName())) {
                    throw new RuntimeException("Unable to match operator: " + fields[1]);
                }
                if (fields.length < 4) {
                    throw new RuntimeException("Operator missing values: " + fields[1]);
                }
                operator.setAcceptCount(Integer.parseInt(fields[2]));
                operator.setRejectCount(Integer.parseInt(fields[3]));
                if (operator instanceof CoercableMCMCOperator) {
                    if (fields.length != 5) {
                        throw new RuntimeException("Coercable operator missing parameter: " + fields[1]);
                    }
                    ((CoercableMCMCOperator)operator).setCoercableParameter(Double.parseDouble(fields[4]));
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

            line = in.readLine();
            fields = line.split("\t");
            // Read in all (possibly more than one) trees
            while (fields[0].equals("tree")) {

                for (Model model : Model.CONNECTED_MODEL_SET) {
                    if (model instanceof TreeModel && fields[1].equals(model.getModelName())) {
                        line = in.readLine();
                        line = in.readLine();
                        fields = line.split("\t");
                        //read number of nodes
                        int nodeCount = Integer.parseInt(fields[0]);
                        double[] nodeHeights = new double[nodeCount];
                        for (int i = 0; i < nodeCount; i++) {
                            line = in.readLine();
                            fields = line.split("\t");
                            nodeHeights[i] = Double.parseDouble(fields[1]);
                        }

                        //on to reading edge information
                        line = in.readLine();
                        line = in.readLine();
                        line = in.readLine();
                        fields = line.split("\t");

                        int edgeCount = Integer.parseInt(fields[0]);

                        int[] parents = new int[edgeCount];
                        for (int i = 0; i < edgeCount; i++){
                            parents[i] = -1;
                        }
                        for (int i = 0; i < edgeCount; i++) {
                            line = in.readLine();
                            if (line != null) {
                                fields = line.split("\t");
                                parents[Integer.parseInt(fields[0])] = Integer.parseInt(fields[1]);
                            }
                        }

                        //perform magic with the acquired information

                        //adopt the loaded tree structure; this does not yet copy the traits on the branches
                        ((TreeModel) model).beginTreeEdit();
                        ((TreeModel) model).adoptTreeStructure(parents, nodeHeights);
                        ((TreeModel) model).endTreeEdit();

                        expectedTreeModelNames.remove(model.getModelName());

                    }
                }

                line = in.readLine();
                if (line != null) {
                    fields = line.split("\t");
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
        }

        return state;
    }

}
