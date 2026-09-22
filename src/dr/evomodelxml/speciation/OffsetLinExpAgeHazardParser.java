package dr.evomodelxml.speciation;

import dr.evomodel.speciation.agedependent.agehazard.OffsetLinExpAgeHazard;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Frederik M. Andersen
 *
 * XML parser for {@link OffsetLinExpAgeHazard}:
 *     h(a) = h0 + (1 + r*gamma*a) * exp(-gamma*a)
 *
 * h0, r and gamma are Parameters of dimension 1 (shared) or K (one per epoch).
 *
 * Per-epoch example (2 epochs):
 * <offsetLinExpAgeHazard>
 *     <h0><parameter value="0.05 0.05" lower="0.0"/></h0>
 *     <r><parameter value="1.5 1.5" lower="0.0"/></r>
 *     <gamma><parameter value="0.3 0.3" lower="0.0"/></gamma>
 * </offsetLinExpAgeHazard>
 */
public class OffsetLinExpAgeHazardParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "offsetLinExpAgeHazard";
    private static final String H0 = "h0";
    private static final String R = "r";
    private static final String GAMMA = "gamma";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter h0    = (Parameter) xo.getElementFirstChild(H0);
        Parameter r     = (Parameter) xo.getElementFirstChild(R);
        Parameter gamma = (Parameter) xo.getElementFirstChild(GAMMA);

        int h0Dim = h0.getDimension();
        int rDim  = r.getDimension();
        int gDim  = gamma.getDimension();
        int maxDim = Math.max(h0Dim, Math.max(rDim, gDim));
        for (int d : new int[]{h0Dim, rDim, gDim}) {
            if (d != 1 && d != maxDim) {
                throw new XMLParseException(
                        "h0 (dim=" + h0Dim + "), r (dim=" + rDim + ") and gamma (dim=" + gDim
                        + ") must each have dimension 1 or the common per-epoch dimension");
            }
        }
        return new OffsetLinExpAgeHazard(h0, r, gamma);
    }

    public String getParserDescription() {
        return "Offset linear-exponential age hazard h(a) = h0 + (1 + r*gamma*a) * exp(-gamma*a). " +
               "h0, r and gamma may each have dimension 1 (shared) or K (one per epoch).";
    }

    public Class getReturnType() {
        return OffsetLinExpAgeHazard.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(H0, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            new ElementRule(R, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            new ElementRule(GAMMA, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
    };
}
