/*
 * ContinuousTreeToKML.java
 *
 * Copyright (C) 2002-2010 BEAST Development Team
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

package dr.app.phylogeography.tools;

import dr.app.util.Arguments;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.io.ImportException;
import jebl.evolution.trees.RootedTree;

import java.io.*;

/**
 * @author Philippe Lemey
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class ContinuousTreeToKML {
    public static void main(String[] args) {


        String inputTreeFile = args[0];
        RootedTree tree = null;
        
        try {
            TreeImporter importer = new NexusImporter(new FileReader(inputTreeFile));
            tree = (RootedTree)importer.importNextTree();
        } catch (ImportException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        double height = 1000000;
        if (args.length > 1) {
            height = Double.parseDouble(args[1]);
        }

        System.out.println("plot height " + height);

        double date = 2000;
        if (args.length > 2) {
            date = Double.parseDouble(args[2]);
        }

        String coordinateName = "location";
        if (args.length > 3) {
            coordinateName = args[3];
        }

        boolean makeTreeSlices = false;
        double[] sliceTimes = null;
        double treeSliceBranchWidth = 3;
        boolean showBranchAtMidPoint = false; // shows complete branch for slice if time is more recent than the branch's midpoint
        if (args.length > 4) {
            try{
                sliceTimes = BetterDiscreteTreeToKML.parseVariableLengthDoubleArray(args[4]);
            } catch (Arguments.ArgumentException ae){
                System.err.println("error reading slice heights");
                ae.printStackTrace();
                return;
            }
//            for (int a = 0; a < sliceTimes.length; a++){
//                sliceTimes[a] = date - sliceTimes[a];
//            }
            makeTreeSlices = true;
        }


        ContinuousKML exporter = new ContinuousKML(tree, args[0], height, date, coordinateName);

        try {
            BufferedWriter out1 = new BufferedWriter(new FileWriter(args[0]+".kml"));
            StringBuffer buffer = new StringBuffer();

            //we write the general tree stuff, but when making slices we do not include everything in the buffer compilation
            exporter.writeTreeToKML();

            if (makeTreeSlices) {
                for (int i = 0; i < sliceTimes.length; i++) {
//                    System.out.println(sliceTimes[i]);
                    exporter.writeTreeToKML(sliceTimes[i], treeSliceBranchWidth, showBranchAtMidPoint);
                }
            }

            exporter.compileBuffer(buffer,makeTreeSlices);
            out1.write(buffer.toString());
            out1.close();
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

    }
}
