/*
 * HiddenLinkageModelParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.MetagenomeData;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.xml.*;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageModelParser extends AbstractXMLObjectParser {

    public static final String LINKAGE_GROUP_COUNT = "linkageGroupCount";

    public static final String NAME = "HiddenLinkageModel";

    
	public String getParserDescription() {
		return "A model to describe missing information about linkage among several reads from a metagenome";
	}


	public Class getReturnType() {
		return HiddenLinkageModel.class;
	}


	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String linkageGroupCount = xo.getAttribute(LINKAGE_GROUP_COUNT, xo.getId());
        MetagenomeData data = (MetagenomeData)xo.getChild(MetagenomeData.class);

        int tc = Integer.parseInt(linkageGroupCount);
        return new HiddenLinkageModel(tc, data);
	}


	public String getParserName() {
		return NAME;
	}

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(LINKAGE_GROUP_COUNT, "The number of hidden lineages", true),
            new ElementRule(MetagenomeData.class)
    };
}
