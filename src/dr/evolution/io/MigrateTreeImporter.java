/*
 * MigrateTreeImporter.java
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

package dr.evolution.io;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;

/**
 * @author Alexei Drummond
 */
public class MigrateTreeImporter extends NexusImporter {

    public static final String POP = "pop";
    public static final String TO_POP = "toPop";
    public static final String FROM_POP = "fromPop";

    public MigrateTreeImporter(Reader reader) {
        super(reader);
    }

    public Tree[] importTrees(TaxonList taxonList) throws IOException, ImportException {
        Tree[] trees = super.importTrees(taxonList);

        for (Tree tree : trees) {

            // convert toPops to pops
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);

                Object toPop = tree.getNodeAttribute(node, TO_POP);
                Object pop = tree.getNodeAttribute(node, POP);

                if (toPop != null && pop != null && !toPop.equals(pop)) {

                    String nodeName = node.getNumber() + "";
                    if (tree.isExternal(node)) {
                        nodeName = tree.getTaxonId(node.getNumber());
                    }
                    if (tree.isRoot(node)) nodeName = "root";

                    throw new RuntimeException(TO_POP + " = " + toPop + ", " + POP + " = " + pop + " in node " + nodeName);
                }

                if (pop == null && toPop != null) {
                    ((MutableTree) tree).setNodeAttribute(node, POP, toPop);
                }
            }
            // convert fromPops to pops
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);

                if (!tree.isRoot(node)) {

                    NodeRef parent = tree.getParent(node);
                    Object fromPop = tree.getNodeAttribute(node, FROM_POP);
                    Object pop = tree.getNodeAttribute(parent, POP);

                    if (fromPop != null && pop != null && !fromPop.equals(pop)) {
                        throw new RuntimeException(FROM_POP + " = " + fromPop + ", " + POP + " = " + pop);
                    }

                    if (pop == null && fromPop != null) {
                        ((MutableTree) tree).setNodeAttribute(parent, POP, fromPop);
                    }
                }
            }

            // fill gaps with parent pops
            fillInternalGaps(tree, tree.getRoot());
            fillExternalGaps(tree);
        }

        return trees;
    }

    private void fillExternalGaps(Tree tree) {

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            if (tree.getNodeAttribute(node, POP) == null) {
                ((MutableTree) tree).setNodeAttribute(node, POP,
                        tree.getNodeAttribute(tree.getParent(node), POP));
            }
        }
    }

    private Object fillInternalGaps(Tree tree, NodeRef node) {

        if (!tree.isExternal(node)) {
            Object left = fillInternalGaps(tree, tree.getChild(node, 0));
            Object right = fillInternalGaps(tree, tree.getChild(node, 1));

            if (tree.getNodeAttribute(node, POP) == null) {

                if (left == null && right == null) {
                    throw new RuntimeException("left and right are both null for node " + node.getNumber());
                }

                if (left == null) left = right;
                if (right == null) right = left;

                if (left.equals(right)) {
                    ((MutableTree) tree).setNodeAttribute(node, POP, left);
                    //System.out.println("Setting pop to " + left + " in node " + node.getNumber());
                } else {
                    throw new RuntimeException(left + "!=" + right + " in children of node " + node.getNumber());
                }
            }
        }
        return tree.getNodeAttribute(node, POP);
    }

    /**
     * Reads a branch in. This could be a node or a tip (calls readNode or readTip
     * accordingly). It then reads the branch length and SimpleNode that will
     * point at the new node or tip.
     */
    FlexibleNode readBranch(HashMap<String, Taxon> translationList) throws IOException, ImportException {
        double length = 0.0;
        FlexibleNode branch;

        clearLastMetaComment();

        if (nextCharacter() == '(') {
            // is an internal node
            branch = readInternalNode(translationList);

        } else {
            // is an external node
            branch = readExternalNode(translationList);
        }

        if (getLastDelimiter() != ':' && getLastDelimiter() != ',' && getLastDelimiter() != ')') {
            String label = readToken(",():;");
            if (label.length() > 0) {
                branch.setAttribute("label", label);
            }
        }

        if (getLastDelimiter() == ':') {
            length = readDouble(" ,():;");

            if (getLastMetaComment() != null) {
                parseMigrationString(getLastMetaComment(), branch);
                clearLastMetaComment();
            }

        }

        branch.setLength(length);

        return branch;
    }

    /**
     * Reads a node in. This could be a polytomy. Calls readBranch on each branch
     * in the node.
     */
    FlexibleNode readInternalNode(HashMap<String, Taxon> translationList) throws IOException, ImportException {
        FlexibleNode node = new FlexibleNode();

        // read the opening '('
        readCharacter();

        // read the first child

        FlexibleNode firstChild = readBranch(translationList);

        node.addChild(firstChild);

        // an internal node must have at least 2 children
        if (getLastDelimiter() != ',') {
            throw new BadFormatException("Missing ',' in tree in TREES block");
        }

        // read subsequent children
        do {
            node.addChild(readBranch(translationList));

        } while (getLastDelimiter() == ',');

        // should have had a closing ')'
        if (getLastDelimiter() != ')') {
            throw new BadFormatException("Missing closing ')' in tree in TREES block");
        }

        readToken(":(),;");

        if (getLastMetaComment() != null) {
            parseMigrationString(getLastMetaComment(), node);
            clearLastMetaComment();
        }

        // find the next delimiter
        return node;
    }

    /**
     * Reads an external node in.
     */
    FlexibleNode readExternalNode(HashMap<String, Taxon> translationList) throws ImportException, IOException {
        FlexibleNode node = new FlexibleNode();

        String label = readToken(":(),;");

        Taxon taxon;

        if (translationList.size() > 0) {
            taxon = translationList.get(label);

            if (taxon == null) {
                // taxon not found in taxon list...
                throw new UnknownTaxonException("Taxon in tree, '" + label + "' is unknown");
            }
        } else {
            taxon = new Taxon(label);
        }

        if (getLastMetaComment() != null) {
            parseMigrationString(getLastMetaComment(), node);
            clearLastMetaComment();
        }

        node.setTaxon(taxon);

        int pop = Integer.parseInt(label.split("\\.")[0]);

        node.setAttribute(POP, (pop - 1));

        return node;
    }

    private void parseMigrationString(String migrationString, FlexibleNode branch) {

        //System.out.println("Parse migration, string='" + migrationString + "'");

        String[] migrations = migrationString.split(";");
        for (String migration : migrations) {
            parseMigration(migration, branch);
        }
    }

    private void parseMigration(String migration, FlexibleNode branch) {

        String[] parts = migration.split(" ");
        if (!parts[0].equals("M")) throw new RuntimeException(parts[0] + " should be M");
        int from = Integer.parseInt(parts[1]);
        parts = parts[2].split(":");
        int to = Integer.parseInt(parts[0]);
        double time = Double.parseDouble(parts[1]);

        if (branch.getAttribute(TO_POP) == null) branch.setAttribute(TO_POP, to);

        branch.setAttribute(FROM_POP, from);
    }


}

