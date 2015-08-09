/*
 * HiddenLinkageTreeLoggerParser.java
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

import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.HiddenLinkageTreeLogger;
import dr.evomodel.tree.TreeLogger;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageTreeLoggerParser extends TreeLoggerParser {

    private XMLSyntaxRule[] rules;

	public HiddenLinkageTreeLoggerParser(){
		rules = new XMLSyntaxRule[super.getSyntaxRules().length + 1];
		System.arraycopy(super.getSyntaxRules(), 0, rules, 0, rules.length-1);
		rules[rules.length-1] = new ElementRule(HiddenLinkageModel.class);
	}
    public String getParserName() {
        return "logHiddenLinkageTree";
    }

    public String getParserDescription() {
        return "Logs a tree with hidden linkage among metagenomic reads to a file";
    }

    public Class getReturnType() {
        return HiddenLinkageTreeLogger.class;
    }
    
    public XMLSyntaxRule[] getSyntaxRules() {
        return this.rules;
    }
    
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	parseXMLParameters(xo);
    	HiddenLinkageModel hlm = (HiddenLinkageModel)xo.getChild(HiddenLinkageModel.class);

    	HiddenLinkageTreeLogger logger = new HiddenLinkageTreeLogger(hlm, tree, branchRates,
                treeAttributeProviders, treeTraitProviders,
                formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition/*,
                normaliseMeanRateTo*/);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;
    }
}
