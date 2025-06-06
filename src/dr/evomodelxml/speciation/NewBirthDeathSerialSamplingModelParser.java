/*
 * NewBirthDeathSerialSamplingModelParser.java
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
import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Arrays;
import java.util.logging.Logger;

public class NewBirthDeathSerialSamplingModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_SERIAL_MODEL = "newBirthDeathSerialSampling";
    public static final String BIRTH_RATE = "birthRate";
    public static final String DEATH_RATE = "deathRate";
    public static final String SAMPLING_RATE = "samplingRate";
    public static final String TREATMENT_PROBABILITY = "treatmentProbability";
    public static final String SAMPLING_PROBABILITY = "samplingProbability";
    public static final String ORIGIN = "origin";
    public static final String CONDITION = "conditionOnSurvival";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String R0 = "R0";
    public static final String D = "D";
    public static final String S = "S";
    public static final String GRADIENT_FLAG = "gradientFlag";
    public static final String GRIDS = "grids";

    public String getParserName() {
        return BIRTH_DEATH_SERIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
        Parameter originParameter = (Parameter) xo.getChild(ORIGIN).getChild(Parameter.class);
        Parameter numGridPoints = (Parameter) xo.getChild(NUM_GRID_POINTS).getChild(Parameter.class);
        Parameter cutoff = (Parameter) xo.getChild(CUT_OFF).getChild(Parameter.class);
        Parameter grids = xo.hasChildNamed(GRIDS) ? (Parameter) xo.getChild(GRIDS).getChild(Parameter.class) : null;

        Parameter lambda = null;
        Parameter mu = null;
        Parameter psi = null;
        Parameter r = null;
        Parameter rho = null;

        Parameter R0Parameter = null;
        Parameter DParameter = null;
        Parameter SParameter = null;

        boolean[] gradientFlags = new boolean[5];
        Arrays.fill(gradientFlags, true);

        if (xo.hasChildNamed(BIRTH_RATE)) {

            lambda = (Parameter) xo.getChild(BIRTH_RATE).getChild(Parameter.class);
            mu = (Parameter) xo.getChild(DEATH_RATE).getChild(Parameter.class);
            psi = (Parameter) xo.getChild(SAMPLING_RATE).getChild(Parameter.class);

            if (xo.getChild(BIRTH_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[0] = xo.getChild(BIRTH_RATE).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(DEATH_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[1] = xo.getChild(DEATH_RATE).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(SAMPLING_RATE).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[2] = xo.getChild(SAMPLING_RATE).getAttribute(GRADIENT_FLAG, true);
            }
        } else {

            R0Parameter = (Parameter) xo.getChild(R0).getChild(Parameter.class);
            DParameter = (Parameter) xo.getChild(D).getChild(Parameter.class);
            SParameter = (Parameter) xo.getChild(S).getChild(Parameter.class);


            if (xo.getChild(R0).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[0] = xo.getChild(R0).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(D).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[1] = xo.getChild(D).getAttribute(GRADIENT_FLAG, true);
            }
            if (xo.getChild(S).hasAttribute(GRADIENT_FLAG)) {
                gradientFlags[2] = xo.getChild(S).getAttribute(GRADIENT_FLAG, true);
        }
        }

        r = (Parameter) xo.getChild(TREATMENT_PROBABILITY).getChild(Parameter.class);
        rho = (Parameter) xo.getChild(SAMPLING_PROBABILITY).getChild(Parameter.class);

        if (xo.getChild(TREATMENT_PROBABILITY).hasAttribute(GRADIENT_FLAG)) {
            gradientFlags[3] = xo.getChild(TREATMENT_PROBABILITY).getAttribute(GRADIENT_FLAG, true);
        }
        if (xo.getChild(SAMPLING_PROBABILITY).hasAttribute(GRADIENT_FLAG)) {
            gradientFlags[4] = xo.getChild(SAMPLING_PROBABILITY).getAttribute(GRADIENT_FLAG, true);
        }

        boolean condition = xo.getAttribute(CONDITION, false);

        NewBirthDeathSerialSamplingModel model;

        if (R0Parameter != null) {
            model = NewBirthDeathSerialSamplingModel.createWithCompoundParameters(
                    R0Parameter, DParameter, SParameter,
                    r, rho, originParameter,
                    condition, (int)numGridPoints.getParameterValue(0), 
                    cutoff.getParameterValue(0), units);
        } else {
            model = new NewBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, originParameter, 
                    condition, (int)(numGridPoints.getParameterValue(0)), 
                    cutoff.getParameterValue(0), units);
        }


        model.setupGradientFlags(gradientFlags);

        model.setupTimeline(grids != null ? grids.getParameterValues() : null);

        Logger.getLogger("dr.evomodel").info("Using birth-death serial sampling model");
        Logger.getLogger("dr.evomodel").info("\tCondition on root: " + condition);
        Logger.getLogger("dr.evomodel").info("\tGrid size: " + (int)numGridPoints.getParameterValue(0));
        Logger.getLogger("dr.evomodel").info("\tGrid end: " + cutoff.getParameterValue(0));
        Logger.getLogger("dr.evomodel").info("\tParameter mode: " + (R0Parameter != null ? "compound" : "raw"));

        return model;
    }

    public String getParserDescription() {
        return "A birth-death model with serial sampling through time and sampled ancestors.";
    }

    public Class getReturnType() {
        return NewBirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(CONDITION, true),
            new XORRule(
                    new ElementRule(BIRTH_RATE,
                            new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class),
                                AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                            }),
                    new ElementRule(R0,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
            ),
            new ElementRule(DEATH_RATE,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }, true),
            new ElementRule(SAMPLING_RATE,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }, true),
            new ElementRule(D,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(S,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(TREATMENT_PROBABILITY,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }),
            new ElementRule(SAMPLING_PROBABILITY,
                    new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                        AttributeRule.newBooleanRule(GRADIENT_FLAG, true)
                    }),
            new ElementRule(ORIGIN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(NUM_GRID_POINTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(CUT_OFF,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GRIDS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}