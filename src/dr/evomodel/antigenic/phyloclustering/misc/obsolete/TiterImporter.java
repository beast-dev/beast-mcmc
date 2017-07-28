package dr.evomodel.antigenic.phyloclustering.misc.obsolete;
/*
 * NewickImporter.java
 *
 * Copyright (C) 2002-2010 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

//package dr.evolution.io;

import dr.evolution.io.Importer;
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

/**
 * Class for importing Newick tree file format
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: NewickImporter.java,v 1.20 2005/12/07 11:25:35 rambaut Exp $
 */
public class TiterImporter extends Importer  {
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
    public TiterImporter(Reader reader) {
        super(reader);
    }

    /**
     * @param treeString a string containing a tree in newick format
     */
    public TiterImporter(String treeString) {
        this(new StringReader(treeString));
    }

    private String virusStrain[];
    private String transformed_titer[];
    private String titer[];

    public void readTiter(int numNodes) throws IOException, ImportException{
    	
    	virusStrain = new String[numNodes];
    	transformed_titer = new String[numNodes];
    	titer = new String[numNodes];

    	readToken();
    	readToken();
    	readToken();
    	readToken();
    	//System.out.println(readToken());
    	//System.out.println(readToken());
    	//System.out.println(readToken());
    	//System.out.println(readToken());

    	for(int i=0; i < (numNodes); i++){
    		//System.out.println(readInteger()); // we just assume 1 to numNodes now
    		readInteger();
    		virusStrain[i] = readToken();
    		//System.out.println(virusStrain[i]);
 		
    		transformed_titer[i] = readToken();
    		//System.out.println(transformed_titer[i]);
  		
    		try{
    		//titer[i] = read() ;
    			titer[i] = readLine();
    			//System.out.println(titer[i]);
    		} catch (EOFException e) {
            //throw new RuntimeException(ite.getMessage());
    			System.out.println("err");
    		}
    		

    	}
   	
    }
    
    public String getVirusStrain(int i) {
		return virusStrain[i];
	}

	public String getTransformed_titer(int i) {
		return transformed_titer[i];
	}

	public String getTiter(int i) {
		return titer[i];
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

            return new FlexibleTree(root, false, true);

        } catch (EOFException e) {
            throw new ImportException("incomplete tree");
        }
    }

    /**
     * importTrees.
     */
    public Tree[] importTrees(TaxonList taxonList) throws IOException, ImportException {
        boolean done = false;
        ArrayList<FlexibleTree> array = new ArrayList<FlexibleTree>();

        do {

            try {

                skipUntil("(");
                unreadCharacter('(');

                FlexibleNode root = readInternalNode(taxonList);
                FlexibleTree tree = new FlexibleTree(root, false, true);
                array.add(tree);

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

        Tree[] trees = new Tree[array.size()];
        array.toArray(trees);

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

        // an internal node must have at least 2 children
        if (getLastDelimiter() != ',') {
            throw new BadFormatException("Expecting ',' in tree, but got '" + (char) getLastDelimiter() + "'");
        }

        // read subsequent children
        do {
            node.addChild(readBranch(taxonList));

        } while (getLastDelimiter() == ',');

        // should have had a closing ')'
        if (getLastDelimiter() != ')') {
            throw new BadFormatException("Missing closing ')' in tree");
        }

        // If there is a label before the colon, store it:
        try {
            String label = readToken(",():;");
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
