/*
 * GlmCovariateImportanceParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.GlmCovariateImportance;
import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class GlmCovariateImportanceParser extends AbstractXMLObjectParser {
    public final static boolean DEBUG = true;

    private static final String PARSER_NAME = "glmCovariateImportance";

    public String getParserName() {
        return PARSER_NAME;
    }

    @SuppressWarnings("deprecation")
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
        OldGLMSubstitutionModel substitutionModel = (OldGLMSubstitutionModel) xo.getChild(OldGLMSubstitutionModel.class);

        return new GlmCovariateImportance(likelihood, substitutionModel);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @SuppressWarnings("deprecation")
    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule("fileName", true),
            AttributeRule.newStringRule("parameterColumn", true),
            AttributeRule.newIntegerRule("burnIn", true),
            new XORRule(
                    new XMLSyntaxRule[] {
                            new ElementRule(TreeDataLikelihood.class),
                            new ElementRule(BeagleTreeLikelihood.class),
                            new ElementRule(CompoundLikelihood.class),
                    }),
            new ElementRule(OldGLMSubstitutionModel.class),
    };

    public String getParserDescription() {
        return "Calculates model deviance for each fixed effect in a phylogeographic GLM";
    }

    public Class getReturnType() {
        return Loggable.class;
    }
}

