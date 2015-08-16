/*
 * HiddenLinkageTreeLogger.java
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

package dr.evomodel.tree;

import java.text.NumberFormat;
import java.util.Set;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.loggers.LogFormatter;

/**
 * @author Aaron Darling (koadman)
 * This class logs hidden linkage trees
 * It modifies a tree to add nodes for metagenomic reads according to a hidden linkage model
 */
public class HiddenLinkageTreeLogger extends TreeLogger {

	HiddenLinkageModel hlm;
	Tree originalTree;
    private LogUpon condition = null;
	
    public HiddenLinkageTreeLogger(HiddenLinkageModel hlm, Tree tree, BranchRates branchRates,
            TreeAttributeProvider[] treeAttributeProviders,
            TreeTraitProvider[] treeTraitProviders,
            LogFormatter formatter, int logEvery, boolean nexusFormat,
            boolean sortTranslationTable, boolean mapNames, NumberFormat format,
            TreeLogger.LogUpon condition) {
    	super(processTree(tree, hlm), branchRates, treeAttributeProviders, treeTraitProviders,formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition);

    	this.originalTree = tree;
    	this.condition = condition;
    	this.hlm = hlm;
    }

    public void log(long state) {
        final boolean doIt = condition != null ? condition.logNow(state) :
            (logEvery < 0 || ((state % logEvery) == 0));
        if(!doIt)
        	return;
        setTree(processTree(originalTree, hlm));
    	super.log(state);
    }
    /**
     * Create a tree of the same topology as the input tree, but with 
     * the linkage groups replaced by their constituent reads
     */
    protected static SimpleTree processTree(Tree tree, HiddenLinkageModel hlm)
    {
    	TaxonList reads = hlm.getData().getReadsTaxa();
    	TaxonList reference = hlm.getData().getReferenceTaxa();
    	
    	// allocate space
    	int nodeCount = tree.getTaxonCount() + reads.getTaxonCount();
    	nodeCount = 2*nodeCount - 1;
    	SimpleNode[] nodes = new SimpleNode[nodeCount];
    	for(int i=0; i<nodes.length; i++){
    		nodes[i] = new SimpleNode();
    		nodes[i].setNumber(i);
    	}
    	SimpleNode root = null;

    	// copy the tree structure
    	for(int i=0; i<tree.getNodeCount(); i++){
    		NodeRef n = tree.getNode(i);
    		for(int cI=0; cI<tree.getChildCount(n); cI++)
    		{
	    		NodeRef c1 = tree.getChild(n, cI);
	    		nodes[n.getNumber()].addChild(nodes[c1.getNumber()]);
    		}
    		nodes[n.getNumber()].setHeight(tree.getNodeHeight(n));
    		nodes[n.getNumber()].setRate(tree.getNodeRate(n));
    		nodes[n.getNumber()].setTaxon(tree.getNodeTaxon(n));
    	}
    	root = nodes[tree.getRoot().getNumber()];
    	
    	// now replace linkage groups with their constituent reads
    	// first free up anything in the range of read leaf nodes
    	int nextFree=tree.getNodeCount();
    	int readI=0;
    	for(int i=reference.getTaxonCount(); i<reference.getTaxonCount()+reads.getTaxonCount(); i++)
    	{
    		SimpleNode tmp = nodes[nextFree];
    		nodes[nextFree] = nodes[i];
    		nodes[nextFree].setNumber(nextFree);
    		nodes[i] = tmp;
    		nodes[i].setNumber(i);
    		nodes[i].setTaxon(reads.getTaxon(readI));
    		readI++;
    		nextFree++;
    	}

    	// now find all linkage group nodes.
    	// if a linkage group has one read, then swap in the read's node
    	// if a linkage group has no reads, delete it and the parent.
    	// if the linkage group has many reads, build a ladder
    	for(int i=0; i<nodes.length; i++)
    	{
    		SimpleNode n = nodes[i];
    		if(n.getTaxon()==null)
    			continue;
    		if(reads.getTaxonIndex(n.getTaxon())>=0 || 
    			reference.getTaxonIndex(n.getTaxon())>= 0)
    			continue;	// not a linkage group

    		int gid = hlm.getTaxonIndex(n.getTaxon()) - reference.getTaxonCount();
    		if(gid<0){
    			System.err.println("big trouble, little china");
    		}
    		Set<Taxon> group = hlm.getGroup(gid);
    		if(group.size()==0){
    			// remove the group completely
    			SimpleNode parent = n.getParent();
    			parent.removeChild(n);
    			if(parent.getChildCount()==1){
        			SimpleNode grandparent = parent.getParent();
    				SimpleNode child = parent.getChild(0);
    				parent.removeChild(child);
    				if(grandparent==null)
    				{
    					root = child;	// parent is root!  other child should become root.
    				}else{
	    				grandparent.removeChild(parent);
	    				grandparent.addChild(child);
    				}
    			}
    		}else if(group.size()==1){
    			// swap the group with the constituent read
    			Taxon[] tax = new Taxon[group.size()];
    			tax = (Taxon[])group.toArray(tax);
    			int rI = getTaxonNode(tax[0], nodes);
    			SimpleNode parent = n.getParent();
    			parent.removeChild(n);
    			parent.addChild(nodes[rI]);
    		}else{
    			// create a star tree with the reads
    			Taxon[] tax = new Taxon[group.size()];
    			tax = (Taxon[])group.toArray(tax);
    			SimpleNode parent = n.getParent();
    			parent.removeChild(n);
    			parent.addChild(nodes[nextFree]);
    			int tI=0;
    			for(; tI < tax.length-2; tI++)
    			{
    				int rI = getTaxonNode(tax[tI], nodes);
    				nodes[nextFree].addChild(nodes[rI]);
    				nodes[nextFree].addChild(nodes[nextFree+1]);
    				nextFree++;
    			}
				int rI = getTaxonNode(tax[tI], nodes);
    			nodes[nextFree].addChild(nodes[rI]);
				int rJ = getTaxonNode(tax[tI+1], nodes);
    			nodes[nextFree].addChild(nodes[rJ]);
    			nextFree++;
    		}
    	}
    	
    	SimpleTree st = new SimpleTree(root);
    	return st;
    }
    
    private static int getTaxonNode(Taxon t, SimpleNode[] nodes){
		int rI=0;
		for(; rI<nodes.length; rI++){
			if(nodes[rI].getTaxon()==t)
				break;
		}
		return rI;
    }
}
