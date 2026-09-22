package dr.evomodelxml.speciation;

import dr.evomodel.speciation.agedependent.agehazard.LinExpAgeHazard;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Frederik M. Andersen
 *
 * XML parser for {@link LinExpAgeHazard}.
 *
 * r and gamma are Parameters of dimension 1 (shared across all epochs)
 * or dimension K (one value per epoch). Mixed dimensions (r dim=1, gamma dim=K
 * or vice versa) are also accepted.
 *
 * Single-epoch example:
 * <linExpAgeHazard>
 *     <r><parameter value="0.5" lower="0.0"/></r>
 *     <gamma><parameter value="0.3" lower="0.0"/></gamma>
 * </linExpAgeHazard>
 *
 * Per-epoch example (2 epochs):
 * <linExpAgeHazard>
 *     <r><parameter value="0.4 0.7" lower="0.0"/></r>
 *     <gamma><parameter value="0.25 0.1" lower="0.0"/></gamma>
 * </linExpAgeHazard>
 */
public class LinExpAgeHazardParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "linExpAgeHazard";
    private static final String R = "r";
    private static final String GAMMA = "gamma";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter r     = (Parameter) xo.getElementFirstChild(R);
        Parameter gamma = (Parameter) xo.getElementFirstChild(GAMMA);

        int rDim = r.getDimension();
        int gDim = gamma.getDimension();
        if (rDim != 1 && gDim != 1 && rDim != gDim) {
            throw new XMLParseException(
                    "r (dim=" + rDim + ") and gamma (dim=" + gDim + ") must either "
                    + "both have the same dimension or one must have dimension 1");
        }

        return new LinExpAgeHazard(r, gamma);
    }

    public String getParserDescription() {
        return "Linear-exponential age hazard h(a) = (1 + r*gamma*a) * exp(-gamma*a). " +
               "r and gamma may each have dimension 1 (shared) or K (one per epoch).";
    }

    public Class getReturnType() {
        return LinExpAgeHazard.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(R, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(GAMMA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };
}
