/*
 * AlloppNode.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.alloppnet.speciation;



import java.util.Stack;

import dr.evolution.tree.SimpleNode;
import dr.evolution.util.Taxon;

import jebl.util.FixedBitSet;


/**
 * An AlloppNode is an interface implemented by DipHistNode in AlloppDiploidHistory, 
 * and MulLabNode in AlloppMulLabTree.
 * The ABC AlloppNode.Abstract contains some common functionality.
 * 
 * @author Graham Jones
 *         Date: 20/03/2012
 */


public interface AlloppNode {
	int nofChildren();
	AlloppNode getChild(int ch);
	AlloppNode getAnc();
	Taxon getTaxon();
	double getHeight();
	FixedBitSet getUnion();
	
	void setChild(int ch, AlloppNode newchild);
	void setAnc(AlloppNode anc);
	void setTaxon(String name);
	void setHeight(double height);
	void setUnion(FixedBitSet union);
	void addChildren(AlloppNode c0, AlloppNode c1);

    String asText(int indentlen);
	
	void fillinUnionsInSubtree(int unionsize);
	AlloppNode nodeOfUnionInSubtree(FixedBitSet x);


    public abstract class Abstract implements AlloppNode {
		 

		public void fillinUnionsInSubtree(int unionsize) {
			if (nofChildren() > 0) {
				getChild(0).fillinUnionsInSubtree(unionsize);
				getChild(1).fillinUnionsInSubtree(unionsize);
				FixedBitSet union = new FixedBitSet(unionsize);
				union.union(getChild(0).getUnion());
				union.union(getChild(1).getUnion());
				setUnion(union);
			}
		}	
		
		
		/* 
		 * For constructor of AlloppDiploidHistory, AlloppMulLabTree.
		 * Searches subtree rooted at node for most tipward node 
		 * whose union contains x. If x is known to be a union of one of the nodes,
		 * it finds that node, so acts as a map union -> node
		 */
		public AlloppNode nodeOfUnionInSubtree(FixedBitSet x) {
			if (nofChildren() == 0) {
				return this;
			}
			if (x.setInclusion(getChild(0).getUnion())) {
				return getChild(0).nodeOfUnionInSubtree(x);
			} else if (x.setInclusion(getChild(1).getUnion())) {
				return getChild(1).nodeOfUnionInSubtree(x);
			} else {
				return this;
			}
		}

        public static String subtreeAsText(AlloppNode node, String s, Stack<Integer> x, int depth, String b) {
            Integer[] y = x.toArray(new Integer[x.size()]);
            StringBuffer indent = new StringBuffer();
            for (int i = 0; i < depth; i++) {
                indent.append("  ");
            }
            for (int i = 0; i < y.length; i++) {
                indent.replace(2*y[i], 2*y[i]+1, "|");
            }
            if (b.length() > 0) {
                indent.replace(indent.length()-b.length(), indent.length(), b);
            }
            s += indent;
            s += node.asText(indent.length());
            s += System.getProperty("line.separator");
            String subs = "";
            if (node.nofChildren() > 0) {
                x.push(depth);
                subs += subtreeAsText(node.getChild(0), "", x, depth+1, "-");
                x.pop();
                subs += subtreeAsText(node.getChild(1), "", x, depth+1, "`-");
            }
            return s + subs;
        }

	}


}
