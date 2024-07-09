/*
 * NewickImporter.java
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

package dr.evolution.io;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.io.EOFException;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for importing Newick tree file format
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class NewickImporter extends Importer implements TreeImporter {
    public static final String COMMENT = "comment";

    public class BranchMissingException extends ImportException {
        /**
         *
         */
        private static final long serialVersionUID = 777435104809244693L;

        public BranchMissingException() {
            super();
        }

        public BranchMissingException(String msg) {
            super("Branch missing: " + msg);
            System.err.println(msg);
        }
    }

    /**
     * @param reader A reader to a source containing a tree in Newick format
     */
    public NewickImporter(Reader reader) {
        super(reader);
    }

    /**
     * @param treeString a string containing a tree in newick format
     */
    public NewickImporter(String treeString) {
        this(new StringReader(treeString));
    }

    /**
     * importTree.
     */
    public Tree importTree(TaxonList taxonList) throws IOException, ImportException {
        setCommentDelimiters('[', ']', '\0', '\0', '&');

        try {
            skipUntil("(");
            unreadCharacter('(');

            final FlexibleNode root = readInternalNode(taxonList);
            if (getLastMetaComment() != null) {
                root.setAttribute(COMMENT, getLastMetaComment());
            }
//			if (getLastDelimiter() != ';') {
//				throw new BadFormatException("Expecting ';' after tree");
//			}

            Map<Taxon, Integer> taxonNumberMap = null;
            if (taxonList != null) {
                taxonNumberMap = new HashMap<Taxon, Integer>();
                int i = 0;
                for (Taxon taxon : taxonList) {
                    taxonNumberMap.put(taxon, i);
                    i++;
                }
            }

            return new FlexibleTree(root, false, true, taxonNumberMap);

        } catch (EOFException e) {
            throw new ImportException("incomplete tree");
        }
    }

    /**
     * countTrees.
     * Counts the number of trees in the file without importing them
     */
    public int countTrees() throws IOException {
        boolean done = false;
        int count = 0;

        do {
            try {
                skipUntil("(");
                skipUntil(";");
                count++;

            } catch (EOFException e) {
                done = true;
            }
        } while (!done);

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

                FlexibleNode root = readInternalNode(taxonList);
                FlexibleTree tree = new FlexibleTree(root, false, true);
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
            skipUntil("(");
            unreadCharacter('(');
        } catch (EOFException e) {
            lastTree = null;
            return false;
        }

        return true;
    }

    private Tree lastTree = null;

    /**
     * import the next tree.
     * return the tree or null if no more trees are available
     */
    public Tree importNextTree() throws IOException, ImportException {
        FlexibleTree tree = null;

        try {
            skipUntil("(");
            unreadCharacter('(');

            FlexibleNode root = readInternalNode(lastTree);

            tree = new FlexibleTree(root, false, true);

        } catch (EOFException e) {
            //
        }

        lastTree = tree;

        return tree;
    }

    /**
     * Reads a branch in. This could be a node or a tip (calls readNode or readTip
     * accordingly). It then reads the branch length and SimpleNode that will
     * point at the new node or tip.
     */
    private FlexibleNode readBranch(TaxonList taxonList) throws IOException, ImportException {
        double length = 0.0;
        FlexibleNode branch;

        if (nextCharacter() == '(') {
            // is an internal node
            branch = readInternalNode(taxonList);

        } else {
            // is an external node
            branch = readExternalNode(taxonList);
        }

        final String comment = getLastMetaComment();
        if (comment != null) {
            branch.setAttribute(COMMENT, comment);
            clearLastMetaComment();
        }

        if (getLastDelimiter() == ':') {
            length = readDouble(",():;");
        }

        branch.setLength(length);

        return branch;
    }

    /**
     * Reads a node in. This could be a polytomy. Calls readBranch on each branch
     * in the node.
     */
    private FlexibleNode readInternalNode(TaxonList taxonList) throws IOException, ImportException {
        FlexibleNode node = new FlexibleNode();

        // read the opening '('
        final char ch = readCharacter();
        assert ch == '(';

        // read the first child
        node.addChild(readBranch(taxonList));

        if (getLastDelimiter() != ',') {
            java.util.logging.Logger.getLogger("dr.evolution.io").warning("Internal node only has a single child.");
        }


        // read subsequent children
        while (getLastDelimiter() == ',') {
            node.addChild(readBranch(taxonList));
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
            if (label.length() > 0) {
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
    private FlexibleNode readExternalNode(TaxonList taxonList) throws IOException, ImportException {
        FlexibleNode node = new FlexibleNode();

        String label = readToken(":(),;");

        Taxon taxon;

        if (taxonList != null) {
            // if a taxon list is given then the taxon must be in it...
            int index = taxonList.getTaxonIndex(label);
            if (index != -1) {
                taxon = taxonList.getTaxon(index);
            } else {
                throw new UnknownTaxonException("Taxon in tree, '" + label + "' is unknown");
            }
        } else {
            // No taxon list given so create new taxa
            taxon = new Taxon(label);
        }

        node.setTaxon(taxon);
        return node;
    }

}
