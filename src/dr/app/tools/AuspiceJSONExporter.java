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

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.io.TreeExporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

import java.awt.*;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AuspiceJSONExporter implements TreeExporter {

    private final PrintStream out;
    private NumberFormat formatter = null;


    public AuspiceJSONExporter(PrintStream out) {
        this.out = out;
    }

    /**
     * @param trees the array of trees to export
     */
    public void exportTrees(Tree[] trees) {
    }

    /**
     * Export a tree with all its attributes.
     *
     * @param tree the tree to export.
     */
    public void exportTree(Tree tree) {
        StringBuilder builder = new StringBuilder();

        double x = 0;
        double[] y = { 0 };

        appendTree(tree, tree.getRoot(), x, y, builder, 0);

        out.print(builder.toString());
    }


    private double appendTree(Tree tree, NodeRef node, double x, double[] y, StringBuilder builder, int indent) {
        builder.append("{\n");

        double x1 = x + tree.getBranchLength(node);
        double y1 = 0.0;

        if (tree.isExternal(node)) {
            appendIndent(builder, indent + 1);
            appendTaxonName(tree, node, builder);
            builder.append(",\n");
            appendAttributes(tree, node, builder, indent + 1);
            y1 = y[0];
            y[0] += 1.0;
        } else {
            appendAttribute("clade", node.getNumber(), builder, indent + 1, false);
            appendIndent(builder, indent + 1);
            builder.append("\"children\": [\n");
            final int last = tree.getChildCount(node) - 1;
            y1 = 0.0;
            for (int i = 0; i < tree.getChildCount(node); i++) {
                appendIndent(builder, indent + 2);
                y1 += appendTree(tree, tree.getChild(node, i), x1, y, builder, indent + 2);
                builder.append(i == last ? "\n" : ",\n");
            }
            y1 /= tree.getChildCount(node);

            appendIndent(builder, indent + 1);
            builder.append("],\n");
        }

//        appendIndent(builder, indent + 1);
//        builder.append("\"height\": ").append(tree.getHeight(node));

        appendAttribute("xvalue", x1, builder, indent + 1, false);
        appendAttribute("yvalue", y1, builder, indent + 1, true);

        appendIndent(builder, indent);
        builder.append("}");

        return y1;
    }

    private StringBuilder appendTaxonName(Tree tree, NodeRef node, StringBuilder builder) {
        String name = tree.getNodeTaxon(node).getId();
        if (name.contains("\"")) {
            name = name.replace("\"", "\\\"");
        }
        return builder.append("\"name\": \"").append(name).append("\"");
    }

    public static double roundDouble(double value, int decimalPlace) {
        double power_of_ten = 1;
        while (decimalPlace-- > 0)
            power_of_ten *= 10.0;
        return Math.round(value * power_of_ten) / power_of_ten;
    }

    /*
	  "num_date": 2015.17803,
	  "xvalue": 0.00229,
	  "country": "SLE",
	  "region": "WesternUrban",
	  "lab": "IG",
	  "strain": "3004_C1_KT5787",
	  "date": "2015-03-07",
	  "freq": {
		"global": "undefined"
	  },
	  "nuc_muts": "",
	  "clade": 17,
	  "yvalue": 15
	*/
    private StringBuilder appendAttributes(Tree tree, NodeRef node, StringBuilder builder, int indent) {
        String id = tree.getNodeTaxon(node).getId();
        String[] fields = id.split("\\|");

        appendAttribute("strain", fields[2], builder, indent, false);
        appendAttribute("lab", fields[3], builder, indent, false);
        appendAttribute("country", fields[4], builder, indent, false);
        appendAttribute("region", fields[5], builder, indent, false);
        appendAttribute("date", fields[8], builder, indent, false);
        return builder;
    }

    private StringBuilder appendAttribute(String key, Object value, StringBuilder builder, int indent, boolean isLast) {
        appendIndent(builder, indent);

        builder.append("\"").append(key).append("\": ");
        appendAttributeValue(value, builder);

        if (!isLast) {
            builder.append(",");
        }
        builder.append("\n");

        return builder;
    }

    private StringBuilder appendAttributeValue(Object value, StringBuilder builder) {
        if (value instanceof Object[]) {
            builder.append("[");
            Object[] elements = ((Object[])value);

            if (elements.length > 0) {
                appendAttributeValue(elements[0], builder);
                for (int i = 1; i < elements.length; i++) {
                    builder.append(",");
                    appendAttributeValue(elements[i], builder);
                }
            }
            return builder.append("]");
        }

        if (value instanceof Color) {
            return builder.append("#").append(((Color)value).getRGB());
        }

        if (value instanceof String) {
            if (((String) value).contains("+")) {
                return builder.append("\"").append(((String) value).split("\\+")[0]).append("\"");
            }
            return builder.append("\"").append(value).append("\"");
        }

        return builder.append(value);
    }

    private void appendIndent(StringBuilder builder, int indent) {
        for (int i = 0; i < indent; i++) {
            builder.append('\t');
        }
    }


}// END: class
