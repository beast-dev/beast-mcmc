/*
 * MascotGradientParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public final class MascotGradientParser extends AbstractXMLObjectParser {

    @Override
    public String getParserName() {
        return MascotGradient.MASCOT_GRADIENT;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MascotLikelihood likelihood = (MascotLikelihood) xo.getChild(MascotLikelihood.class);
        MascotGradient gradient = new MascotGradient(likelihood);
        gradient.setNumericGradientStepSize(NumericGradientStepSizeProvider.parseStepSizeRatio(xo));
        return gradient;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MascotLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return "Gradient of a MASCOT likelihood with respect to its flat log-parameter vector.";
    }

    @Override
    public Class getReturnType() {
        return MascotGradient.class;
    }
}
