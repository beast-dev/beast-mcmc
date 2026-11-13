
package dr.evomodelxml.treedatalikelihood;

import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.DiscretePreOrderReport;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;
import dr.xml.*;

/**
 * Reportable XML helper for dumping all cached discrete pre-order partials.
 */
public class DiscretePreOrderReportParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "discretePreOrderReport";
    private static final String TOLERANCE = "tolerance";
    private static final String EIGEN_BASIS = "eigenBasis";
    private static final String TYPE = "type";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        if (likelihood == null) {
            throw new XMLParseException("Expected a treeDataLikelihood child in " + PARSER_NAME);
        }
        double tolerance = xo.getAttribute(TOLERANCE, 0.0);

        String typeString = xo.getAttribute(TYPE, DiscretePartialsType.TOP.getMeaning());
        DiscretePartialsType type = DiscretePartialsType.parse(typeString);

        return new DiscretePreOrderReport(likelihood, type, tolerance);
    }

    @Override
    public String getParserDescription() {
        return "Reports branch-start and branch-end pre-order partials for a discrete TreeDataLikelihood.";
    }

    @Override
    public Class getReturnType() {
        return DiscretePreOrderReport.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newDoubleRule(TOLERANCE, true),
                AttributeRule.newStringRule(TYPE, true),
        };
    }

    public static final XMLObjectParser PARSER = new DiscretePreOrderReportParser();
}
