/*
 * NexusExporter.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.io.TreeExporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.Iterator;
import java.util.Map;
import java.util.HashMap;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: NexusExporter.java,v 1.5 2006/09/08 11:34:53 rambaut Exp $
 */
public class NexusExporter implements TreeExporter {

    public NexusExporter(PrintStream out) {
        this.out = out;
    }

    public void exportTrees(Tree[] trees) {
        writeNexusHeader(trees[0]);
        out.println("\t\t;");
        for (int i = 0; i < trees.length; i++) {
            writeNexusTree(trees[i], i);
        }
        out.println("End;");
    }

    public void exportTree(Tree tree) {
        writeNexusHeader(tree);
        out.println("\t\t;");
        writeNexusTree(tree, 1);
        out.println("End;");
    }

    private void writeNexusTree(Tree tree, int i) {
        out.print("tree TREE" + i  + "  = [&R] ");
        writeNode( tree, tree.getRoot(), true, null);
        out.println(";");
    }

    private Map<String, Integer> writeNexusHeader(Tree tree) {
        int taxonCount = tree.getTaxonCount();
        out.println("#NEXUS");
        out.println();
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + taxonCount + ";");
        out.println("\tTaxlabels");
        for (int i = 0; i < taxonCount; i++) {
            out.println("\t\t" + tree.getTaxon(i).getId());
        }
        out.println("\t\t;");
        out.println("End;");
        out.println("");
        out.println("Begin trees;");

        // This is needed if the trees use numerical taxon labels
        out.println("\tTranslate");
        Map<String, Integer> idMap = new HashMap<String, Integer>();
        for (int i = 0; i < taxonCount; i++) {
            int k = i + 1;
            idMap.put(tree.getTaxonId(i),k);
            if (k < taxonCount) {
                out.println("\t\t" + k + " " + tree.getTaxonId(i) + ",");
            } else {
                out.println("\t\t" + k + " " + tree.getTaxonId(i));
            }
        }
        return idMap;
    }

    private void writeNode(Tree tree, NodeRef node, boolean attributes, Map<String, Integer> idMap) {
        if (tree.isExternal(node)) {
            int k = node.getNumber() + 1;
            if (idMap != null) k = idMap.get(tree.getTaxonId(k-1));

            out.print(k);
        } else {
            out.print("(");
            writeNode(tree, tree.getChild(node, 0), attributes, idMap);
            for (int i = 1; i < tree.getChildCount(node); i++) {
                out.print(",");
                writeNode(tree, tree.getChild(node, i), attributes, idMap);
            }
            out.print(")");
        }

        Iterator iter = tree.getNodeAttributeNames(node);
        if (iter != null) {
            boolean first = true;
            while (iter.hasNext()) {
                if (first) {
                    out.print("[&");
                    first = false;
                } else {
                    out.print(",");
                }
                String name = (String)iter.next();
                out.print(name + "=");
                Object value = tree.getNodeAttribute(node, name);
                printValue(value);
            }
            out.print("]");
        }

        if (!tree.isRoot(node)) {
            out.print(":");

            double length = tree.getBranchLength(node);
            out.print(Double.toString(length));
        }
    }

    private void printValue(Object value) {
        if (value instanceof Object[]) {
            out.print("{");
            Object[] values = (Object[])value;
            for (int i = 0; i < values.length; i++) {
                if (i > 0) {
                    out.print(",");
                }
                printValue(values[i]);
            }
            out.print("}");
        } else if (value instanceof String) {
            out.print("\"" + value.toString() + "\"");
        } else {
            out.print(value.toString());
        }
    }

    private PrintStream out;

}
