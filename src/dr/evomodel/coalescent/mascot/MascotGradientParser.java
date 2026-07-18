/*
 * MascotGradientParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public final class MascotGradientParser extends AbstractXMLObjectParser {

    public static final String PART = "part";
    public static final String MIGRATION = "migration";
    public static final String POPULATION_SIZE = "popSizes";
    public static final String CLOCK_RATE = "clockRate";

    @Override
    public String getParserName() {
        return MascotGradient.MASCOT_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MascotLikelihood likelihood = (MascotLikelihood) xo.getChild(MascotLikelihood.class);
        String partString = xo.getStringAttribute(PART);
        MascotGradient.Part part;
        if (partString.equals(MIGRATION)) {
            part = MascotGradient.Part.MIGRATION;
        } else if (partString.equals(POPULATION_SIZE)) {
            part = MascotGradient.Part.POPULATION_SIZE;
        } else if (partString.equals(CLOCK_RATE)) {
            part = MascotGradient.Part.CLOCK_RATE;
        } else {
            throw new XMLParseException("mascotGradient part must be '" + MIGRATION + "', '" +
                    POPULATION_SIZE + "', or '" + CLOCK_RATE + "', found '" + partString + "'");
        }
        MascotGradient gradient = new MascotGradient(likelihood, part);
        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(PART),
            new ElementRule(MascotLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return "Gradient of a MASCOT likelihood with respect to its migration-rate, population-size, " +
                "or clock-rate parameter (part=\"migration\", part=\"popSizes\", or part=\"clockRate\").";
    }

    @Override
    public Class getReturnType() {
        return MascotGradient.class;
    }
}
