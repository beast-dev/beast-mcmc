/*
 * NexusExporter.java
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

import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.io.TreeExporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: NexusExporter.java,v 1.5 2006/09/08 11:34:53 rambaut Exp $
 */
public class NexusExporter implements TreeExporter {

    private final PrintStream out;
    private NumberFormat formatter = null;
    private String treePrefix = DEFAULT_TREE_PREFIX;
    private boolean sorted = false;
    private AttributeType writeAttributesAs = AttributeType.NODE_ATTRIBUTES;
    public static final String DEFAULT_TREE_PREFIX = "TREE";
    public static final String SPECIAL_CHARACTERS_REGEX = ".*[\\s\\.;,\"\'].*";
    
    public enum AttributeType {
        NODE_ATTRIBUTES,
        BRANCH_ATTRIBUTES
    }

//    /*
//     * this constructor will result in npe if calling methods other than exportAlignment()
//     * */
//    public NexusExporter() {
//        this.out = null;
//    }
    
    public NexusExporter(PrintStream out) {
        this.out = out;
        this.writeAttributesAs = AttributeType.NODE_ATTRIBUTES;
    }

    public NexusExporter(PrintStream out, AttributeType writeAttributesAs) {
        this.out = out;
        this.writeAttributesAs = writeAttributesAs;
    }

    /**
     * Sets the name to use for each tree (will be suffixed by tree number)
     *
     * @param treePrefix
     */
    public void setTreePrefix(String treePrefix) {
        this.treePrefix = treePrefix;
    }

    /**
     * Sets the number format to use for outputting branch lengths
     *
     * @param format
     */
    public void setNumberFormat(NumberFormat format) {
        formatter = format;
    }

    /**
     * @param sorted true if you wish the translation table to be alphabetically sorted.
     */
    public void setSortedTranslationTable(boolean sorted) {
        this.sorted = sorted;
    }

    /**
     * @param trees      the array of trees to export
     * @param attributes true if the nodes should be annotated with their attributes
     * @param treeNames  Names of the trees
     */
    public void exportTrees(Tree[] trees, boolean attributes, String[] treeNames) {
        if(!(treeNames==null) && trees.length != treeNames.length) {
            throw new RuntimeException("Number of trees and number of tree names is not the same");
        }
        Map<String, Integer> idMap = writeNexusHeader(trees[0]);
        out.println("\t\t;");
        for (int i = 0; i < trees.length; i++) {
            if(treeNames==null) {
                writeNexusTree(trees[i], treePrefix + i, attributes, idMap);
            }
            else {
                writeNexusTree(trees[i], treeNames[i], attributes, idMap);
            }
        }
        out.println("End;");
    }

    public void exportTrees(Tree[] trees, boolean attributes) {
        exportTrees(trees, attributes, null);
    }


    public void exportTrees(Tree[] trees) {
        exportTrees(trees, true, null);
    }

    /**
     * Export a tree with all its attributes.
     *
     * @param tree the tree to export.
     */
    public void exportTree(Tree tree) {
        Map<String, Integer> idMap = writeNexusHeader(tree);
        out.println("\t\t;");
        writeNexusTree(tree, treePrefix + 1, true, idMap);
        out.println("End;");
    }

    public void writeNexusTree(Tree tree, String s, boolean attributes, Map<String, Integer> idMap) {
        // PAUP marks rooted trees thou
        String treeAttributes = "[&R] ";

        // Place tree level attributes in tree comment
        StringBuilder treeComment = null;
        {
            Iterator<String> iter = tree.getAttributeNames();
            if (iter != null) {
                while (iter.hasNext()) {
                    final String name = iter.next();
                    final String value = tree.getAttribute(name).toString();

                    if( name.equals("weight") ) {
                        treeAttributes = treeAttributes + "[&W " + value + " ] ";
                    }
                    else {
                        if( treeComment == null ) {
                            treeComment = new StringBuilder(" [&");
                        } else if( treeComment.length() > 2 ) {
                            treeComment.append(", ");
                        }

                        treeComment.append(name).append("=").append(value);
                    }
                }
                if( treeComment != null ) {
                    treeComment.append("]");
                }
            }
        }

        out.print("tree " + s + ((treeComment != null) ? treeComment.toString() : "")
                + " = " + treeAttributes);
        
        writeNode(tree, tree.getRoot(), attributes, idMap);
        out.println(";");
    }

    public Map<String, Integer> writeNexusHeader(Tree tree) {
        int taxonCount = tree.getTaxonCount();
        List<String> names = new ArrayList<String>();

        for (int i = 0; i < tree.getTaxonCount(); i++) {
            names.add(tree.getTaxonId(i));
        }

        if (sorted) Collections.sort(names);

        out.println("#NEXUS");
        out.println();
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + taxonCount + ";");
        out.println("\tTaxlabels");
        for (String name : names) {
            if (name.matches(SPECIAL_CHARACTERS_REGEX)) {
                name = "'" + name + "'";
            }
            out.println("\t\t" + name);
        }
        out.println("\t\t;");
        out.println("End;");
        out.println("");
        out.println("Begin trees;");

        // This is needed if the trees use numerical taxon labels
        out.println("\tTranslate");
        Map<String, Integer> idMap = new HashMap<String, Integer>();

        int k = 1;
        for (String name : names) {
            idMap.put(name, k);
            if (name.matches(SPECIAL_CHARACTERS_REGEX)) {
                name = "'" + name + "'";
            }
            if (k < names.size()) {
                out.println("\t\t" + k + " " + name + ",");
            } else {
                out.println("\t\t" + k + " " + name);
            }
            k += 1;
        }
        return idMap;
    }

    private void writeNode(Tree tree, NodeRef node, boolean attributes, Map<String, Integer> idMap) {
        if (tree.isExternal(node)) {
            int k = node.getNumber() + 1;
            if (idMap != null) k = idMap.get(tree.getTaxonId(k - 1));

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

        if (writeAttributesAs == AttributeType.BRANCH_ATTRIBUTES && !tree.isRoot(node)) {
            out.print(":");
        }

        if (attributes) {
            Iterator<?> iter = tree.getNodeAttributeNames(node);
            if (iter != null) {
                boolean first = true;
                while (iter.hasNext()) {
                    if (first) {
                        out.print("[&");
                        first = false;
                    } else {
                        out.print(",");
                    }
                    String name = (String) iter.next();
                    out.print(name + "=");
                    Object value = tree.getNodeAttribute(node, name);
                    printValue(value);
                }
                out.print("]");
            }
        }

        if (writeAttributesAs == AttributeType.NODE_ATTRIBUTES && !tree.isRoot(node)) {
            out.print(":");
        }

        if (!tree.isRoot(node)) {
            double length = tree.getBranchLength(node);
            if (formatter != null) {
                out.print(formatter.format(length));
            } else {
                out.print(length);
            }
        }
    }

    private void printValue(Object value) {
        if (value instanceof Object[]) {
            out.print("{");
            Object[] values = (Object[]) value;
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

	public String exportAlignment(Alignment alignment) throws IOException, IllegalArgumentException {

		StringBuffer buffer = new StringBuffer();
		
		DataType dataType = null;
		int seqLength = 0;
		
		for (int i = 0; i < alignment.getSequenceCount(); i++) {

			Sequence sequence = alignment.getSequence(i);

			if (sequence.getLength() > seqLength) {
				seqLength = sequence.getLength();
			}

			if (dataType == null) {
				dataType = sequence.getDataType();
			} else if (dataType != sequence.getDataType()) {
				throw new RuntimeException(
						"Sequences must have the same data type.");
			}// END: dataType check

		}// END: sequences loop

		buffer.append("#NEXUS\n");
		buffer.append("begin data;\n");
		buffer.append("\tdimensions" + " " + "ntax=" + alignment.getTaxonCount() + " " + "nchar=" + seqLength + ";\n");
		buffer.append("\tformat datatype=" + dataType.getDescription()
				+ " missing=" + DataType.UNKNOWN_CHARACTER + " gap="
				+ DataType.GAP_CHARACTER + ";\n");
		buffer.append("\tmatrix\n");
		
		int maxRowLength = seqLength;
		for (int n = 0; n < Math.ceil((double) seqLength / maxRowLength); n++) {
			for (int i = 0; i < alignment.getSequenceCount(); i++) {

				Sequence sequence = alignment.getSequence(i);

				StringBuilder builder = new StringBuilder("\t");

				appendTaxonName(sequence.getTaxon(), builder);

				String sequenceString = sequence.getSequenceString();

				builder.append("\t").append(
						sequenceString.subSequence(
								n * maxRowLength,
								Math.min((n + 1) * maxRowLength,
										sequenceString.length())));
				int shortBy = Math.min(Math.min(n * maxRowLength, seqLength)
						- sequence.getLength(), maxRowLength);
				
				if (shortBy > 0) {
					for (int j = 0; j < shortBy; j++) {
						builder.append(DataType.GAP_CHARACTER);
					}
				}
				
				buffer.append(builder + "\n");
			}// END: sequences loop
		}
		buffer.append(";\nend;");
		
		return buffer.toString();
	}// END: exportAlignment

    /**
     * name suitable for printing - quotes if necessary
     * @param taxon
     * @param builder
     * @return
     */
    private StringBuilder appendTaxonName(Taxon taxon, StringBuilder builder) {
        
    	String name = taxon.getId();
        if (!name.matches(SPECIAL_CHARACTERS_REGEX)) {
            // JEBL way of quoting the quote character
            name = name.replace("\'", "\'\'");
            builder.append("\'").append(name).append("\'");
            return builder;
        }
        return builder.append(name);
    }
	
}// END: class
