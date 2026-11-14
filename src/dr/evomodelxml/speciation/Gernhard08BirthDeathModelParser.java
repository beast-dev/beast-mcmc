/*
 * BirthDeathModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.Gernhard08BirthDeathModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Joseph Heled
 */
public class Gernhard08BirthDeathModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_MODEL = "Gernhard08BirthDeathModel";
    public static final String BIRTHDIFF_RATE = "birthMinusDeathRate";
    public static final String RELATIVE_DEATH_RATE = "relativeDeathRate";
    public static final String SAMPLE_PROB = "sampleProbability";
    public static final String TREE_TYPE = "type";
    public static final String CONDITIONAL_ON_ROOT = "conditionalOnRoot";

    public static final String BIRTH_DEATH = "birthDeath";
    public static final String MEAN_GROWTH_RATE_PARAM_NAME = BIRTH_DEATH + ".meanGrowthRate";
    public static final String RELATIVE_DEATH_RATE_PARAM_NAME = BIRTH_DEATH + "." + RELATIVE_DEATH_RATE;

    public String getParserName() {
        return BIRTH_DEATH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        final String s = xo.getAttribute(TREE_TYPE, Gernhard08BirthDeathModel.TreeType.UNSCALED.toString());
        final Gernhard08BirthDeathModel.TreeType treeType = Gernhard08BirthDeathModel.TreeType.valueOf(s);
        final boolean conditonalOnRoot =  xo.getAttribute(CONDITIONAL_ON_ROOT, false);

        final Parameter birthParameter = (Parameter) xo.getElementFirstChild(BIRTHDIFF_RATE);
        final Parameter deathParameter = (Parameter) xo.getElementFirstChild(RELATIVE_DEATH_RATE);
        final Parameter sampleProbability = xo.hasChildNamed(SAMPLE_PROB) ?
                (Parameter) xo.getElementFirstChild(SAMPLE_PROB) : null;

        Logger.getLogger("dr.evomodel").info(xo.hasChildNamed(SAMPLE_PROB) ? getCitationRHO() : getCitation());

        final String modelName = xo.getId();

        return new Gernhard08BirthDeathModel(modelName, birthParameter, deathParameter, sampleProbability,
                treeType, units, conditonalOnRoot);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
    public static String getCitationRHO() {
        return "Stadler, T; On incomplete sampling under birth-death models and connections to the sampling-based coalescent;\n" +
               "JOURNAL OF THEORETICAL BIOLOGY (2009) 261:58-66";
    }

    public static String getCitation() {
        return "Using birth-death model on tree: Gernhard T (2008) J Theor Biol, Volume 253, Issue 4, Pages 769-778 In press";
    }

    public String getParserDescription() {
        return "Gernhard (2008) model of speciation (equation at bottom of page 19 of draft).";
    }

    public Class getReturnType() {
        return Gernhard08BirthDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TREE_TYPE, true),
            AttributeRule.newBooleanRule(CONDITIONAL_ON_ROOT, true),
            new ElementRule(BIRTHDIFF_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RELATIVE_DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SAMPLE_PROB, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}