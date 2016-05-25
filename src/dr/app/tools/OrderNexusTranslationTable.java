/*
 * OrderNexusTranslationTable.java
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

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;

import java.io.IOException;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.text.NumberFormat;

/**
 * A simple routine to alphabetically order the taxa translation table
 * in a nexus trees file.
 *
 * @author Alexei Drummond
 */
public class OrderNexusTranslationTable {

    public static void main(String[] args) throws Importer.ImportException, IOException {

        NexusImporter importer = new NexusImporter(new FileReader(args[0]));

        Tree[] trees = importer.importTrees(null);
        System.out.println("Read " + trees.length + " trees from " + args[0]);

        String newFileName = args[0] + ".new";

        PrintStream ps = new PrintStream(new FileOutputStream(newFileName));

        NexusExporter exporter = new NexusExporter(ps);

        exporter.setTreePrefix("STATE_");
        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(7);
        exporter.setNumberFormat(format);
        exporter.setSortedTranslationTable(true);

        exporter.exportTrees(trees, false);

        ps.flush();
        ps.close();
        System.out.println("Wrote " + trees.length + " trees to " + newFileName);
    }

}
