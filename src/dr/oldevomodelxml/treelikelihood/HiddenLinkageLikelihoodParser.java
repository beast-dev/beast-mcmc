/*
 * HiddenLinkageLikelihoodParser.java
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

package dr.oldevomodelxml.treelikelihood;

import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.HiddenLinkageLikelihood;
import dr.xml.*;

/**
 * @author Aaron Darling (koadman)
 */
@Deprecated // Switching to BEAGLE
public class HiddenLinkageLikelihoodParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "A likelihood calculator for hidden linkage among metagenomic reads";
	}

	@Override
	public Class getReturnType() {
		return HiddenLinkageLikelihood.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		HiddenLinkageModel hlm = (HiddenLinkageModel) xo.getChild(HiddenLinkageModel.class);
		TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        return new HiddenLinkageLikelihood(hlm, tree);
	}


	public String getParserName() {
		return "HiddenLinkageLikelihood";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(HiddenLinkageModel.class),
            new ElementRule(TreeModel.class),
    };

}
