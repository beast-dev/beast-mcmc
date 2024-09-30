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
import dr.evomodel.speciation.MasBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.TwoParamBirthDeathSerialSamplingModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class NewBirthDeathSerialSamplingModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_SERIAL_MODEL = "newBirthDeathSerialSampling";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String LAMBDA = "birthRate";
    public static final String MU = "deathRate";
    public static final String PSI = "samplingRate";
    public static final String RHO = "samplingProbability";
    public static final String R = "treatmentProbability";
    public static final String ORIGIN = "origin";
    public static final String TREE_TYPE = "type";
    public static final String CONDITION = "conditionOnSurvival";

    public static final String GRADIENT_FLAG = "gradientFlag";

    public static final String GRIDS = "grids";

    public static final String BDSS = "bdss";

    enum ParametersToUse {
        ALL,
        LAM_MU_PSI,
        LAM_PSI
    }

    public String getParserName() {
        return BIRTH_DEATH_SERIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        final Parameter cutoff = xo.hasChildNamed(CUT_OFF) ? (Parameter) xo.getElementFirstChild(CUT_OFF) : new Parameter.Default(Double.POSITIVE_INFINITY);
        final Parameter numGridPoints = xo.hasChildNamed(NUM_GRID_POINTS) ? (Parameter) xo.getElementFirstChild(NUM_GRID_POINTS): new Parameter.Default(1.0);

        final Parameter lambda = (Parameter) xo.getElementFirstChild(LAMBDA);
        final Parameter mu = (Parameter) xo.getElementFirstChild(MU);
        final Parameter psi = (Parameter) xo.getElementFirstChild(PSI);

        final Parameter rho    = xo.hasChildNamed(RHO) ? (Parameter) xo.getElementFirstChild(RHO) : new Parameter.Default(0.0);
        final Parameter r      = xo.hasChildNamed(R) ? (Parameter) xo.getElementFirstChild(R) : new Parameter.Default(1.0);

        final Parameter origin = (Parameter) xo.getElementFirstChild(ORIGIN);

        Boolean condition = xo.getAttribute(CONDITION, false);

        final boolean[] gradientFlags = new boolean[5];

        gradientFlags[0] = xo.getChild(LAMBDA).getAttribute(GRADIENT_FLAG, true).booleanValue();
        gradientFlags[1] = xo.getChild(MU).getAttribute(GRADIENT_FLAG, true).booleanValue();
        gradientFlags[2] = xo.getChild(PSI).getAttribute(GRADIENT_FLAG, true).booleanValue();
        gradientFlags[3] = xo.getChild(RHO).getAttribute(GRADIENT_FLAG, true).booleanValue();
        gradientFlags[4] = xo.getChild(R).getAttribute(GRADIENT_FLAG, true).booleanValue();

        ParametersToUse parametersToUse;

        if (gradientFlags[0] && gradientFlags[1] && gradientFlags[2] && (!gradientFlags[3]) && (!gradientFlags[4])) {
            // three true, two false
            parametersToUse = ParametersToUse.LAM_MU_PSI;
        } else if (gradientFlags[0] && (!gradientFlags[1]) && gradientFlags[2] && (!gradientFlags[3]) && (!gradientFlags[4])) {
            // two true, three false
            parametersToUse = ParametersToUse.LAM_PSI;
        } else {
            parametersToUse = ParametersToUse.ALL;
        }

        final Parameter grids = xo.hasChildNamed(GRIDS) ? (Parameter) xo.getElementFirstChild(GRIDS): null;


        String citeThisModel;
        if ( r.getParameterValue(0) < Double.MIN_VALUE ) {
            citeThisModel = getCitationFBD();
        } else if ( rho.getParameterValue(0) > Double.MIN_VALUE ) {
            citeThisModel = getCitationSABDP();
        } else {
            citeThisModel = getCitationR0();
        }

        Logger.getLogger("dr.evomodel").info(citeThisModel);

//        NewBirthDeathSerialSamplingModel model = MAS_TEST ?
//                new MasBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, origin, condition, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units) :
//                new NewBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, origin, condition, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units);

        NewBirthDeathSerialSamplingModel model = null;

        switch (parametersToUse) {
            case ALL:
                model = new NewBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, origin, condition, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units);
                break;
            case LAM_MU_PSI:
                model = new MasBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, origin, condition, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units);
                break;
            case LAM_PSI:
                model = new TwoParamBirthDeathSerialSamplingModel(lambda, mu, psi, r, rho, origin, condition, (int)(numGridPoints.getParameterValue(0)), cutoff.getParameterValue(0), units);
                break;
        }

        // model.setupGradientFlags(gradientFlags);
        model.setupTimeline(grids != null ? grids.getParameterValues(): null);
        return model;
    }

    private static final boolean MAS_TEST = true;

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static String getCitationFBD() {
        return "Stadler, T; Sampling-through-time in birth–death trees, " +
                "J. Theor. Biol. (2010) doi: 10.1016/j.jtbi.2010.09.010";
    }

    public static String getCitationSABDP() {
        return "Gavryushkina et al, Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration, " +
                "PLoS Comp. Biol., doi: 10.1371/journal.pcbi.1003919, 2014";
    }

    public static String getCitationR0() {
        return "Stadler et al, Estimating the basic reproductive number from viral sequence data, " +
                "Mol. Biol. Evol., doi: 10.1093/molbev/msr217, 2012";
    }

    public String getParserDescription() {
        return "A serially-sampled birth-death model with the possibility of treatment and sampling at present.";
    }

    public Class getReturnType() {
        return NewBirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TREE_TYPE, true),
            new ElementRule(ORIGIN, Parameter.class, "The origin of the infection, x0 > tree.rootHeight", false),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PSI, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(R, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(RHO, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(CONDITION, new XMLSyntaxRule[]{new ElementRule(boolean.class)}, true),
            new ElementRule(GRADIENT_FLAG, new XMLSyntaxRule[]{new ElementRule(boolean.class)}, true),
            XMLUnits.SYNTAX_RULES[0]
    };
}