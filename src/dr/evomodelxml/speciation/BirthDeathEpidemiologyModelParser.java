/*
 * BirthDeathEpidemiologyModelParser.java
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
 * @author Andrew Rambaut
 * @author Tanja Stadler
 * @author Alexei Drummond
 * @author Joseph Heled
 */
public class BirthDeathEpidemiologyModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_EPIDEMIOLOGY = "birthDeathEpidemiology";
    public static final String R0 = "R0";
    public static final String RECOVERY_RATE = "recoveryRate";
    public static final String SAMPLING_PROBABILITY = "samplingProbability";
    public static final String ORIGIN = BirthDeathSerialSamplingModelParser.ORIGIN;

    public String getParserName() {
        return BIRTH_DEATH_EPIDEMIOLOGY;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        final Parameter R0Parameter = (Parameter) xo.getElementFirstChild(R0);
        final Parameter recoveryRateParameter = (Parameter) xo.getElementFirstChild(RECOVERY_RATE);
        final Parameter samplingProbabiltyParameter = (Parameter) xo.getElementFirstChild(SAMPLING_PROBABILITY);

        Parameter origin = null;
        if (xo.hasChildNamed(ORIGIN)) {
            origin = (Parameter) xo.getElementFirstChild(ORIGIN);
        }

        Logger.getLogger("dr.evomodel").info("Using epidemiological parameterization of " + getCitationRT());

        return new BirthDeathSerialSamplingModel(modelName, R0Parameter, recoveryRateParameter, samplingProbabiltyParameter, origin, units);
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
        return "Stadler et al (2011) model of epidemiology.";
    }

    public Class getReturnType() {
        return BirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(ORIGIN, Parameter.class, "The origin of the infection, x0 > tree.rootHeight", true),
            new ElementRule(R0, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RECOVERY_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SAMPLING_PROBABILITY, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}