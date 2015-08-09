/*
 * DebugUtils.java
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

package dr.inference.mcmc;

/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.io.*;

public class DebugUtils {

    /**
     * Writes out the current state in a human readable format to help debugging.
     * If it fails, then returns false but does not stop.
     * @param file the file
     * @param state the current state number
     * @return success
     */
    public static boolean writeStateToFile(File file, long state, double lnL) {
        OutputStream fileOut = null;
        try {
            fileOut = new FileOutputStream(file);
            PrintStream out = new PrintStream(fileOut);

            out.print("state\t");
            out.println(state);

            out.print("lnL\t");
            out.println(lnL);

            for (Parameter parameter : Parameter.FULL_PARAMETER_SET) {
                out.print(parameter.getParameterName());
                out.print("\t");
                out.print(parameter.getDimension());
                for (int dim = 0; dim < parameter.getDimension(); dim++) {
                    out.print("\t");
                    out.print(parameter.getParameterValue(dim));
                }
                out.println();
            }

            for (Model model : Model.FULL_MODEL_SET) {
                if (model instanceof TreeModel) {
                    out.print(model.getModelName());
                    out.print("\t");
                    out.println(((TreeModel) model).getNewick());
                }
            }

            out.close();
            fileOut.close();
        } catch (IOException ioe) {
            System.err.println("Unable to write file: " + ioe.getMessage());
            return false;
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
    public static long readStateFromFile(File file, double[] lnL) {
        long state = -1;

        try {
            FileReader fileIn = new FileReader(file);
            BufferedReader in = new BufferedReader(fileIn);

            String line = in.readLine();
            String[] fields = line.split("\t");
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

            for (Parameter parameter : Parameter.FULL_PARAMETER_SET) {
                line = in.readLine();
                fields = line.split("\t");
//                if (!fields[0].equals(parameter.getParameterName())) {
//                    System.err.println("Unable to match state parameter: " + fields[0] + ", expecting " + parameter.getParameterName());
//                }
                int dimension = Integer.parseInt(fields[1]);

                if (dimension != parameter.getDimension()) {
                    System.err.println("Unable to match state parameter dimension: " + dimension + ", expecting " + parameter.getDimension());
                }

                for (int dim = 0; dim < parameter.getDimension(); dim++) {
                    parameter.setParameterValueQuietly(dim, Double.parseDouble(fields[dim + 2]));
                }
            }

            // load the tree models last as we get the node heights from the tree (not the parameters which
            // which may not be associated with the right node
            for (Model model : Model.FULL_MODEL_SET) {
                if (model instanceof TreeModel) {
                    line = in.readLine();
                    fields = line.split("\t");
                    if (!fields[0].equals(model.getModelName())) {
                        throw new RuntimeException("Unable to match state parameter: " + fields[0] + ", expecting " + model.getModelName());
                    }
                    NewickImporter importer = new NewickImporter(fields[1]);
                    Tree tree = importer.importNextTree();
                    ((TreeModel) model).beginTreeEdit();
                    ((TreeModel) model).adoptTreeStructure(tree);
                    ((TreeModel) model).endTreeEdit();
                }
            }

            in.close();
            fileIn.close();
            for (Likelihood likelihood : Likelihood.FULL_LIKELIHOOD_SET) {
                likelihood.makeDirty();
            }
        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
        } catch (Importer.ImportException ie) {
            throw new RuntimeException("Unable to import tree: " + ie.getMessage());
        }

        return state;
    }


}
