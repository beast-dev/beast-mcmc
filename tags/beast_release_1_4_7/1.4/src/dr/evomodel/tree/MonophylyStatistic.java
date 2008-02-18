/*
 * MonophylyStatistic.java
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
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.inference.model.BooleanStatistic;
import dr.xml.*;

import java.util.Set;

/**
 * Performs monophyly test given a taxonList
 *
 * @version $Id: MonophylyStatistic.java,v 1.16 2005/07/11 14:06:25 rambaut Exp $
 *
 * @author Roald Forsberg
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 */
public class MonophylyStatistic extends BooleanStatistic implements TreeStatistic {
       
	public static final String MONOPHYLY_STATISTIC = "monophylyStatistic";
	public static final String MRCA = "mrca";

    public MonophylyStatistic(String name, Tree tree, TaxonList taxa) throws Tree.MissingTaxonException {

		super(name);
		this.tree = tree;
		this.leafSet = Tree.Utils.getLeavesForTaxa(tree, taxa);
	
    }
    
	public void setTree(Tree tree) { this.tree = tree; }
	public Tree getTree() { return tree; }
	
	public int getDimension() { return 1; }
	    
    /** @return boolean result of test. */
    public boolean getBoolean(int dim) {
	
		return Tree.Utils.isMonophyletic(this.tree, this.leafSet);
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	
		public String getParserName() { return MONOPHYLY_STATISTIC; }
	
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}
			Tree tree = (Tree)xo.getChild(Tree.class);
			
			XMLObject cxo = (XMLObject)xo.getChild(MRCA);
			TaxonList taxa = (TaxonList)cxo.getChild(TaxonList.class);
						
			try {
				return new MonophylyStatistic(name, tree, taxa);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
			}
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A statistic that returns true if a given set of taxa are monophyletic for a given tree";
		}

		public Class getReturnType() { return MonophylyStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new StringAttributeRule(NAME, "A name for this statistic for the purpose of logging", true),
			new ElementRule(TreeModel.class),
			new ElementRule(MRCA, new XMLSyntaxRule[] {
				new ElementRule(Taxa.class)
			}),
		};

	};

	private Tree tree = null;
	private Set leafSet = null;
    
}
