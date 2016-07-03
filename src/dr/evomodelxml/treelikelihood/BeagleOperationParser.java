/*
 * BeagleOperationParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treelikelihood;

import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.treelikelihood.BeagleOperationReport;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.io.PrintWriter;

/**
 * @author Marc Suchard
 * @version $Id$
 */
public class BeagleOperationParser extends AbstractXMLObjectParser {


    public static final String OPERATION_REPORT = "beagleOperationReport";
    public static final String BRANCH_FILE_NAME = "branchFileName";
    public static final String OPERATION_FILE_NAME = "operationFileName";

    public String getParserName() {
        return OPERATION_REPORT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SitePatterns patternList = (SitePatterns) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        GammaSiteRateModel substitutionModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        PrintWriter branch = null, operation = null;

        if (xo.hasAttribute(BRANCH_FILE_NAME)) {
            branch = XMLParser.getFilePrintWriter(xo, OPERATION_REPORT, BRANCH_FILE_NAME);
        }
        if (xo.hasAttribute(OPERATION_FILE_NAME)) {
            operation = XMLParser.getFilePrintWriter(xo, OPERATION_REPORT, OPERATION_FILE_NAME);
        }

        return new BeagleOperationReport(treeModel, patternList, rateModel, substitutionModel, alignment, branch, operation);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return BeagleOperationReport.class;
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class, 2, 2),
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class),
            new ElementRule(GammaSiteRateModel.class),
//            new ElementRule(Alignment.class),
            AttributeRule.newStringRule(BRANCH_FILE_NAME, true),
            AttributeRule.newStringRule(OPERATION_FILE_NAME, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
