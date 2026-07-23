/*
 * BEASTTreesImporter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.tools.oldtreeannotator;

import dr.evolution.io.Importer;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.util.*;

/**
 * Class for importing Newick-type tree file format with numbered tips
 * This is intended only for reading BEAST trees files quickly.
 *
 * @author Andrew Rambaut
 */
public class BEASTTreesImporter extends Importer implements TreeImporter {
    public static final String COMMENT = "comment";
    private String treeName = null;

    /**
     * @param reader A reader to a source containing a tree in Newick format
     */
    public BEASTTreesImporter(Reader reader, boolean hasComments) {
        super(reader);
        if (hasComments) {
            setCommentDelimiters('[', ']', '\0', '\0', '&');
        }
    }

    /**
     * importTree.
     */
    public Tree importTree(TaxonList ignored) throws IOException, ImportException {
        throw new UnsupportedOperationException("this method is not available");
    }

    /**
     * countTrees.
     * Counts the number of trees in the file without importing them
     */
    public int countTrees() throws IOException {
        int count = 0;
        String line;

        BufferedReader br = (BufferedReader)getReader();
        do {
            line = br.readLine();
            if (line != null && line.trim().startsWith("tree ")) {
                count++;
            }

        } while (line != null);

        return count;
    }

    /**
     * importTrees.
     */
    public List<Tree> importTrees(TaxonList taxonList) throws IOException, ImportException {
        boolean done = false;
        List<Tree> trees = new ArrayList<>();

        do {

            try {

                skipUntil("(");
                unreadCharacter('(');

                FlexibleNode root = readInternalNode();
                FlexibleTree tree = new FlexibleTree(root);
                trees.add(tree);

                if (taxonList == null) {
                    taxonList = tree;
                }

                if (readCharacter() != ';') {
                    throw new BadFormatException("Expecting ';' after tree");
                }

            } catch (EOFException e) {
                done = true;
            }
        } while (!done);

        return trees;
    }

    /**
     * return whether another tree is available.
     */
    public boolean hasTree() throws IOException, ImportException {
        try {
            skipUntilToken("STATE", "_");
            treeName = readToken(",():;");

            skipUntil("(");
            unreadCharacter('(');
        } catch (EOFException e) {
            return false;
        }

        return true;
    }

    /**
     * import the next tree.
     * return the tree or null if no more trees are available
     */
    public Tree importNextTree() throws IOException, ImportException {
        if (treeName == null) {
            skipUntilToken("STATE", "_");
            treeName = readToken(",():;");
        }

        FlexibleTree tree = null;

        try {
            skipUntil("(");
            unreadCharacter('(');

            FlexibleNode root = readInternalNode();

            tree = new FlexibleTree(root);

        } catch (EOFException e) {
            //
        }

        tree.setId(treeName);
        treeName = null;

        return tree;
    }

    /**
     * Reads a branch in. This could be a node or a tip (calls readNode or readTip
     * accordingly). It then reads the branch length and SimpleNode that will
     * point at the new node or tip.
     */
    private FlexibleNode readBranch() throws IOException, ImportException {
        double length = 0.0;
        FlexibleNode branch;

        if (nextCharacter() == '(') {
            // is an internal node
            branch = readInternalNode();

        } else {
            // is an external node
            branch = readExternalNode();
        }

        final String comment = getLastMetaComment();
        if (comment != null) {
            branch.setAttribute(COMMENT, comment);
            clearLastMetaComment();
        }

        if (getLastDelimiter() == ':') {
            length = readDouble(",():;");
        }

        branch.setAttribute("_length", length);

        return branch;
    }

    /**
     * Reads a node in. This could be a polytomy. Calls readBranch on each branch
     * in the node.
     */
    private FlexibleNode readInternalNode() throws IOException, ImportException {
        FlexibleNode node = new FlexibleNode();

        // read the opening '('
        final char ch = readCharacter();
        assert ch == '(';

        // read the first child
        node.addChild(readBranch());

        if (getLastDelimiter() != ',') {
            java.util.logging.Logger.getLogger("dr.evolution.io").warning("Internal node only has a single child.");
        }


        // read subsequent children
        while (getLastDelimiter() == ',') {
            node.addChild(readBranch());
        }

        // should have had a closing ')'
        if (getLastDelimiter() != ')') {
            throw new BadFormatException("Missing closing ')' in tree");
        }

        // If there is a label before the colon, store it:
        try {
            String label = readToken(",():;");
            if ((char) getLastDelimiter() == ';') {
                unreadCharacter(';');
            }
            if (!label.isEmpty()) {
                node.setAttribute("label", label);
            }
        } catch (IOException ioe) {
            // probably an end of file without a terminal ';'
            // we are going to allow this and return the nodes...
        }

        return node;
    }

    /**
     * Reads an external node in.
     */
    private FlexibleNode readExternalNode() throws IOException {
        FlexibleNode node = new FlexibleNode();
        String label = readToken(":(),;");
        node.setNumber(Integer.parseInt(label) - 1); // node numbers index from 0
        return node;
    }

    public String getTreeName() {
        return treeName;
    }
}
