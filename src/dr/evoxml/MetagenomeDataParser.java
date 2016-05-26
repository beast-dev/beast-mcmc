/*
 * MetagenomeDataParser.java
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

package dr.evoxml;

import dr.evolution.LinkageConstraints;
import dr.evolution.MetagenomeData;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.xml.*;

/**
 * @author Aaron Darling
 */
public class MetagenomeDataParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "Data representing metagenome reads aligned to reference sequences and a reference tree";
	}

	@Override
	public Class getReturnType() {
		return MetagenomeData.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TaxonList taxa = null;
    	taxa = (TaxonList)xo.getChild(TaxonList.class);
        if(taxa==null)
        	taxa = (Tree)xo.getChild(Tree.class);

        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        LinkageConstraints lc = (LinkageConstraints)xo.getChild(LinkageConstraints.class);
        
        boolean fixedReferenceTree = false;
        if (xo.hasAttribute("fixedReferenceTree")) {
        	fixedReferenceTree = xo.getBooleanAttribute("fixedReferenceTree");
        }

        MetagenomeData md = new MetagenomeData(taxa, alignment, lc, fixedReferenceTree); 
        return md;
	}

	
	public String getParserName() {
		return "MetagenomeData";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new XORRule(new ElementRule(Taxa.class),
						new ElementRule(Tree.class)),
            new ElementRule(Alignment.class),
            new ElementRule(LinkageConstraints.class, true),	// optional element
	        AttributeRule.newBooleanRule("fixedReferenceTree", true, "Whether the reference tree should be of fixed topology"),
    };

}
