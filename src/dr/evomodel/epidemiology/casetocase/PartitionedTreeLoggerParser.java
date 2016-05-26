/*
 * PartitionedTreeLoggerParser.java
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

package dr.evomodel.epidemiology.casetocase;

import dr.evomodelxml.tree.TreeLoggerParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * @author Matthew Hall
 */
public class PartitionedTreeLoggerParser extends TreeLoggerParser {

    private XMLSyntaxRule[] rules;

    public PartitionedTreeLoggerParser(){
        rules = new XMLSyntaxRule[super.getSyntaxRules().length + 1];
        System.arraycopy(super.getSyntaxRules(), 0, rules, 0, rules.length-1);
        rules[rules.length-1] = new ElementRule(CaseToCaseTreeLikelihood.class);
    }
    public String getParserName() {
        return "logPartitionedTree";
    }

    public String getParserDescription() {
        return "Logs a partitioned tree (phylogenetic tree and transmission tree)";
    }

    public Class getReturnType() {
        return PartitionedTreeLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return this.rules;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        parseXMLParameters(xo);
        CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)xo.getChild(CaseToCaseTreeLikelihood.class);

        PartitionedTreeLogger logger = new PartitionedTreeLogger(c2cTL, tree, branchRates,
                treeAttributeProviders, treeTraitProviders,
                formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;
    }
}

