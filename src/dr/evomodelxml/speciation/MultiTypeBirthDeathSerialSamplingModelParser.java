/*
 * MultiTypeBirthDeathSerialSamplingModelParser.java
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
import dr.evomodel.speciation.MultiTypeBirthDeathSerialSamplingModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * XML parser for the PROTOTYPE {@link MultiTypeBirthDeathSerialSamplingModel}.
 *
 * Example (single epoch; see {@code numGridPoints}/{@code cutOff}, used exactly as in
 * {@link NewBirthDeathSerialSamplingModelParser}, for the episodic/multi-epoch case):
 * <pre>{@code
 * <multiTypeBirthDeathSerialSampling id="mtbd" types="A B" typeAttribute="state"
 *                                     conditionOnSurvival="true">
 *     <origin><parameter value="6.0"/></origin>
 *     <birthRate><parameter value="2.0 0.1 0.1 1.5"/></birthRate>
 *     <migrationRate><parameter value="0.0 0.2 0.2 0.0"/></migrationRate>
 *     <deathRate><parameter value="0.5 0.5"/></deathRate>
 *     <samplingRate><parameter value="0.3 0.3"/></samplingRate>
 *     <samplingProbability><parameter value="0.5 0.5"/></samplingProbability>
 * </multiTypeBirthDeathSerialSampling>
 * }</pre>
 * birthRate/migrationRate are flattened K*K (row-major, parent-type-major: index i*K+j for
 * birthRate_ij / migrationRate_ij), optionally K*K per epoch (one block per epoch, epochs
 * outermost); deathRate/samplingRate/treatmentProbability/samplingProbability are length K,
 * optionally K per epoch. Unlike {@link NewBirthDeathSerialSamplingModelParser}, a bare
 * scalar is not accepted as shorthand for "the same rate in every type" -- give all K (or
 * K*K) values explicitly.
 *
 * @author Marc A. Suchard
 * @author Claude
 */
public class MultiTypeBirthDeathSerialSamplingModelParser extends AbstractXMLObjectParser {

    public static final String MULTI_TYPE_BIRTH_DEATH_SERIAL_MODEL = "multiTypeBirthDeathSerialSampling";

    public static final String TYPES = "types";
    public static final String TYPE_ATTRIBUTE = "typeAttribute";
    public static final String STEPS_PER_UNIT_TIME = "stepsPerUnitTime";

    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";

    public static final String LAMBDA = "birthRate";
    public static final String MIGRATION = "migrationRate";
    public static final String MU = "deathRate";
    public static final String PSI = "samplingRate";
    public static final String R = "treatmentProbability";
    public static final String RHO = "samplingProbability";
    public static final String ORIGIN = "origin";
    public static final String CONDITION = "conditionOnSurvival";

    public String getParserName() {
        return MULTI_TYPE_BIRTH_DEATH_SERIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        String typesString = xo.getStringAttribute(TYPES);
        String[] typeLabels = typesString.trim().split("[\\s,]+");
        int numTypes = typeLabels.length;
        if (numTypes < 1) {
            throw new XMLParseException("'" + TYPES + "' must list at least one type");
        }

        String typeAttribute = xo.getAttribute(TYPE_ATTRIBUTE, MultiTypeBirthDeathSerialSamplingModel.DEFAULT_TYPE_ATTRIBUTE);
        int stepsPerUnitTime = xo.getAttribute(STEPS_PER_UNIT_TIME, 200);

        final int numGridPoints = xo.hasChildNamed(NUM_GRID_POINTS) ?
                (int) ((Parameter) xo.getElementFirstChild(NUM_GRID_POINTS)).getParameterValue(0) : 1;
        final double cutOff = xo.hasChildNamed(CUT_OFF) ?
                ((Parameter) xo.getElementFirstChild(CUT_OFF)).getParameterValue(0) : Double.POSITIVE_INFINITY;

        final Parameter birthRate = (Parameter) xo.getElementFirstChild(LAMBDA);
        final Parameter migrationRate = (Parameter) xo.getElementFirstChild(MIGRATION);
        final Parameter deathRate = (Parameter) xo.getElementFirstChild(MU);
        final Parameter psi = (Parameter) xo.getElementFirstChild(PSI);

        final Parameter r = xo.hasChildNamed(R) ? (Parameter) xo.getElementFirstChild(R)
                : onesParameter(numTypes);
        final Parameter rho = xo.hasChildNamed(RHO) ? (Parameter) xo.getElementFirstChild(RHO)
                : zerosParameter(numTypes);

        final Parameter origin = (Parameter) xo.getElementFirstChild(ORIGIN);

        boolean condition = xo.getAttribute(CONDITION, true);

        return new MultiTypeBirthDeathSerialSamplingModel(modelName, numTypes, typeLabels,
                birthRate, migrationRate, deathRate, psi, r, rho, origin, condition,
                numGridPoints, cutOff, typeAttribute, stepsPerUnitTime, units);
    }

    private static Parameter onesParameter(int n) {
        double[] v = new double[n];
        java.util.Arrays.fill(v, 1.0);
        return new Parameter.Default(v);
    }

    private static Parameter zerosParameter(int n) {
        return new Parameter.Default(n, 0.0);
    }

    public String getParserDescription() {
        return "A PROTOTYPE episodic multi-type birth-death-sampling model.";
    }

    public Class getReturnType() {
        return MultiTypeBirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TYPES),
            AttributeRule.newStringRule(TYPE_ATTRIBUTE, true),
            AttributeRule.newIntegerRule(STEPS_PER_UNIT_TIME, true),
            AttributeRule.newBooleanRule(CONDITION, true),
            new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(CUT_OFF, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(ORIGIN, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MIGRATION, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PSI, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(R, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(RHO, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}
