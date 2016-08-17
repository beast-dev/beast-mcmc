/*
 * EmpiricalTreeDistributionOperatorParser.java
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

package dr.evomodelxml.operators;

import dr.util.FileHelpers;
import dr.xml.*;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.operators.EmpiricalTreeDistributionOperator;
import dr.inference.operators.MCMCOperator;

import java.io.*;

/**
 * @author Andrew Rambaut
 *         <p/>
 *         Reads a list of trees from a NEXUS file.
 */
public class EmpiricalTreeDistributionOperatorParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return EmpiricalTreeDistributionOperator.EMPIRICAL_TREE_DISTRIBUTION_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        boolean metropolisHastings = false;

        if (xo.hasAttribute(EmpiricalTreeDistributionOperator.METROPOLIS_HASTINGS)) {
            metropolisHastings = xo.getBooleanAttribute(EmpiricalTreeDistributionOperator.METROPOLIS_HASTINGS);
        }

        final EmpiricalTreeDistributionModel treeModel = (EmpiricalTreeDistributionModel) xo.getChild(EmpiricalTreeDistributionModel.class);

        return new EmpiricalTreeDistributionOperator(treeModel, metropolisHastings, weight);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(EmpiricalTreeDistributionOperator.METROPOLIS_HASTINGS, true),
            new ElementRule(EmpiricalTreeDistributionModel.class)
    };

    public String getParserDescription() {
        return "Operator which switches between trees in an empirical distribution.";
    }

    public Class getReturnType() {
        return EmpiricalTreeDistributionOperator.class;
    }
}
