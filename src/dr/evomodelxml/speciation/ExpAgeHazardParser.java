package dr.evomodelxml.speciation;

import dr.evomodel.speciation.agedependent.agehazard.ExpAgeHazard;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Frederik M. Andersen
 *
 * XML parser for {@link ExpAgeHazard}:
 *     h(a) = exp(-gamma*a)
 *
 * gamma is a Parameter of dimension 1 (shared across all epochs) or K (one value
 * per epoch). Positive gamma decays from 1 at age 0; gamma = 0 gives a constant-rate
 * process and negative gamma an increasing hazard.
 *
 * Single-epoch example:
 * <expAgeHazard>
 *     <gamma><parameter value="0.3" lower="0.0"/></gamma>
 * </expAgeHazard>
 *
 * Per-epoch example (2 epochs):
 * <expAgeHazard>
 *     <gamma><parameter value="0.25 0.1" lower="0.0"/></gamma>
 * </expAgeHazard>
 */
public class ExpAgeHazardParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "expAgeHazard";
    private static final String GAMMA = "gamma";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter gamma = (Parameter) xo.getElementFirstChild(GAMMA);
        return new ExpAgeHazard(gamma);
    }

    public String getParserDescription() {
        return "Exponential age hazard h(a) = exp(-gamma*a).";
    }

    public Class getReturnType() {
        return ExpAgeHazard.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules  = {
            new ElementRule(GAMMA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };
}
