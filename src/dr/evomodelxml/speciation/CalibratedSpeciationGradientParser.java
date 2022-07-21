/*
 * CalibratedSpeciationGradientParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.speciation.CalibratedSpeciationGradient;
import dr.evomodel.speciation.CalibratedSpeciationLikelihood;
import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CalibratedSpeciationGradientParser extends AbstractXMLObjectParser {

    public static final String CALIBRATED_SPECIATION_GRADIENT = "calibratedSpeciationGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SpeciationLikelihoodGradient unCalibratedSpeciationGradient = (SpeciationLikelihoodGradient) xo.getChild(SpeciationLikelihoodGradient.class);
        CalibratedSpeciationLikelihood calibratedSpeciationLikelihood = (CalibratedSpeciationLikelihood) xo.getChild(CalibratedSpeciationLikelihood.class);
        return new CalibratedSpeciationGradient(unCalibratedSpeciationGradient, calibratedSpeciationLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SpeciationLikelihoodGradient.class),
            new ElementRule(CalibratedSpeciationLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CalibratedSpeciationGradient.class;
    }

    @Override
    public String getParserName() {
        return CALIBRATED_SPECIATION_GRADIENT;
    }
}
