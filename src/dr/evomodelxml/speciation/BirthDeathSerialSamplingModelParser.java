/*
 * BirthDeathModelParser.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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
    public static final String SAMPLE_PROBABILITY = "sampleProbability";
    public static final String SAMPLED_REMAIN_INFECTIOUS = "sampledRemainInfectiousProb";
    public static final String FINAL_TIME_INTERVAL = "finalTimeInterval";
    public static final String ORIGIN = "origin";
    public static final String TREE_TYPE = "type";
    public static final String BDSS = "bdss";

    public String getParserName() {
        return BIRTH_DEATH_SERIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String modelName = xo.getId();
        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        final Parameter lambda = (Parameter) xo.getElementFirstChild(LAMBDA);

        boolean relativeDeath = xo.hasChildNamed(RELATIVE_MU);
        Parameter mu;
        if (relativeDeath) {
            mu = (Parameter) xo.getElementFirstChild(RELATIVE_MU);
        } else {
            mu = (Parameter) xo.getElementFirstChild(MU);
        }

        final Parameter psi = (Parameter) xo.getElementFirstChild(PSI);
        final Parameter p = (Parameter) xo.getElementFirstChild(SAMPLE_PROBABILITY);

        Parameter origin = null;
        if (xo.hasChildNamed(ORIGIN)) {
            origin = (Parameter) xo.getElementFirstChild(ORIGIN);
        }

        //    for r=1, this is sampledRemainInfectiousProb=0
        //    for r=0, this is sampledRemainInfectiousProb=1
        final Parameter r = xo.hasChildNamed(SAMPLED_REMAIN_INFECTIOUS) ?
                (Parameter) xo.getElementFirstChild(SAMPLED_REMAIN_INFECTIOUS) : new Parameter.Default(1.0);        
        // the r parameter is 1 - sampleRemainsInfectiousProbability
        r.setParameterValueQuietly(0, 1 - r.getParameterValue(0));

        final Parameter finalTimeInterval = xo.hasChildNamed(FINAL_TIME_INTERVAL) ?
                        (Parameter) xo.getElementFirstChild(FINAL_TIME_INTERVAL) : new Parameter.Default(0.0);

        Logger.getLogger("dr.evomodel").info(xo.hasChildNamed(SAMPLED_REMAIN_INFECTIOUS) ? getCitationRT() : getCitationPsiOrg());

        return new BirthDeathSerialSamplingModel(modelName, lambda, mu, psi, p, relativeDeath,
                r, finalTimeInterval, origin, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
    public static String getCitationPsiOrg() {
//        return "Stadler, T; Sampling-through-time in birth-death trees; JOURNAL OF THEORETICAL BIOLOGY (2010) 267:396-404";
        return "Stadler T (2010) J Theor Biol 267, 396-404 [Yule Process with Serial Samples].";
    }
    
    public static String getCitationRT() {
        return "Stadler et al (2011) : Estimating the basic reproductive number from viral sequence data, Submitted.";
    }

    public String getParserDescription() {
        return "Stadler et al (2010; in press) model of speciation.";
    }

    public Class getReturnType() {
        return BirthDeathSerialSamplingModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TREE_TYPE, true),
//            AttributeRule.newDoubleRule(FINAL_TIME_INTERVAL, true),
            new ElementRule(FINAL_TIME_INTERVAL, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(ORIGIN, Parameter.class, "The origin of the infection, x0 > tree.rootHeight", true),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new XORRule(
                    new ElementRule(MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(RELATIVE_MU, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})),
            new ElementRule(PSI, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SAMPLED_REMAIN_INFECTIOUS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(SAMPLE_PROBABILITY, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}