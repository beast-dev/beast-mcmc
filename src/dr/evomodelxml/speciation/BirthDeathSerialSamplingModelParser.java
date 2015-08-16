/*
 * BirthDeathSerialSamplingModelParser.java
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

package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathSerialSamplingModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Joseph Heled
 */
public class BirthDeathSerialSamplingModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_SERIAL_MODEL = "birthDeathSerialSampling";
    public static final String LAMBDA = "birthRate";
    public static final String MU = "deathRate";
    public static final String RELATIVE_MU = "relativeDeathRate";
    public static final String PSI = "psi";
    public static final String SAMPLE_PROBABILITY = "sampleProbability"; // default to fix to 0
    public static final String SAMPLE_BECOMES_NON_INFECTIOUS = "sampleBecomesNonInfectiousProb";
    public static final String R = "r";
//    public static final String FINAL_TIME_INTERVAL = "finalTimeInterval";
    public static final String ORIGIN = "origin";
    public static final String TREE_TYPE = "type";
    public static final String BDSS = "bdss";
    public static final String HAS_FINAL_SAMPLE = "hasFinalSample";

    public String getParserName() {
        return BIRTH_DEATH_SERIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
        
        boolean hasFinalSample = xo.getAttribute(HAS_FINAL_SAMPLE, false);

        final Parameter lambda = (Parameter) xo.getElementFirstChild(LAMBDA);

        boolean relativeDeath = xo.hasChildNamed(RELATIVE_MU);
        Parameter mu;
        if (relativeDeath) {
            mu = (Parameter) xo.getElementFirstChild(RELATIVE_MU);
        } else {
            mu = (Parameter) xo.getElementFirstChild(MU);
        }

        final Parameter psi = (Parameter) xo.getElementFirstChild(PSI);
        //Issue 656: fix p=0
        final Parameter p = xo.hasChildNamed(SAMPLE_PROBABILITY) ?
                        (Parameter) xo.getElementFirstChild(SAMPLE_PROBABILITY) : new Parameter.Default(0.0);

        Parameter origin = null;
        if (xo.hasChildNamed(ORIGIN)) {
            origin = (Parameter) xo.getElementFirstChild(ORIGIN);
        }

        final Parameter r = xo.hasChildNamed(SAMPLE_BECOMES_NON_INFECTIOUS) ?
                (Parameter) xo.getElementFirstChild(SAMPLE_BECOMES_NON_INFECTIOUS) : new Parameter.Default(0.0);
//        r.setParameterValueQuietly(0, 1 - r.getParameterValue(0)); // donot use it, otherwise log is changed improperly

//        final Parameter finalTimeInterval = xo.hasChildNamed(FINAL_TIME_INTERVAL) ?
//                (Parameter) xo.getElementFirstChild(FINAL_TIME_INTERVAL) : new Parameter.Default(0.0);

        Logger.getLogger("dr.evomodel").info(xo.hasChildNamed(SAMPLE_BECOMES_NON_INFECTIOUS) ? getCitationRT() : getCitationPsiOrg());

        return new BirthDeathSerialSamplingModel(modelName, lambda, mu, psi, p, relativeDeath,
                r, hasFinalSample, origin, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public static String getCitationPsiOrg() {
//        return "Stadler, T; Sampling-through-time in birth-death trees; JOURNAL OF THEORETICAL BIOLOGY (2010) 267:396-404";
        return "Stadler T (2010) J Theor Biol 267, 396-404 [Birth-Death with Serial Samples].";
    }

    public static String getCitationRT() {
        return "Stadler et al (2011) : Estimating the basic reproductive number from viral sequence data, " +
                "Mol.Biol.Evol., doi: 10.1093/molbev/msr217, 2011";
    }

    public String getParserDescription() {
        return "Stadler et al (2010) model of speciation.";
    }

    public Class getReturnType() {
        return BirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TREE_TYPE, true),
            AttributeRule.newBooleanRule(HAS_FINAL_SAMPLE, true),
//            new ElementRule(FINAL_TIME_INTERVAL, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(ORIGIN, Parameter.class, "The origin of the infection, x0 > tree.rootHeight", true),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new XORRule(
                    new ElementRule(MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(RELATIVE_MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})),
            new ElementRule(PSI, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SAMPLE_BECOMES_NON_INFECTIOUS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            //Issue 656
//            new ElementRule(SAMPLE_PROBABILITY, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//            XMLUnits.SYNTAX_RULES[0]
    };
}