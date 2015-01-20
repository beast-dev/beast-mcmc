/*
 * TreeMetricStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.inference.model.Statistic;
import dr.xml.*;


/**
 * A statistic that returns the distance between two trees. Currently,
 * this just returns a 0 for identity.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: TreeMetricStatistic.java,v 1.14 2005/07/11 14:06:25 rambaut Exp $
 */
public class TreeMetricStatistic extends Statistic.Abstract implements TreeStatistic {

	public static final String TREE_METRIC_STATISTIC = "treeMetricStatistic";
	public static final String TARGET = "target";
	public static final String REFERENCE = "reference";
	
	public TreeMetricStatistic(String name, Tree target, Tree reference) {

		super(name);
		
		this.target = target;
		this.referenceNewick = Tree.Utils.uniqueNewick(reference, reference.getRoot());

	}

	public void setTree(Tree tree) { this.target = tree; }
	public Tree getTree() { return target; }
	
	public int getDimension() { return 1; }
	
	/** @return value. */
	public double getStatisticValue(int dim) {
		
		if (Tree.Utils.uniqueNewick(target, target.getRoot()).equals(referenceNewick)) {
			return 0.0;
		} else {
			return 1.0;
		}
		
//		return TreeMetrics.getRobinsonFouldsDistance(reference, target);
		
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	
		public String getParserName() { return TREE_METRIC_STATISTIC; }
	
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}
			Tree target = (Tree)xo.getSocketChild(TARGET);
			Tree reference = (Tree)xo.getSocketChild(REFERENCE);
			
			return new TreeMetricStatistic(name, target, reference);
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A statistic that returns the distance between two trees. Currently this is simply a 0 for identity and a 1 for difference.";
		}

		public Class getReturnType() { return TreeMetricStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new StringAttributeRule(NAME, "A name for this statistic primarily for the purposes of logging", true),
			new ElementRule(TARGET, 
				new XMLSyntaxRule[] { new ElementRule(TreeModel.class) }),
			new ElementRule(REFERENCE, 
				new XMLSyntaxRule[] { new ElementRule(Tree.class) }),
		};

	};
	
	private Tree target = null;
	private String referenceNewick = null;
}