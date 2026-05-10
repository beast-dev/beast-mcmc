
package dr.evomodelxml.treedatalikelihood;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.DiscretePreOrderReport;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Reportable XML helper for dumping all cached discrete pre-order partials.
 */
public class DiscretePreOrderReportParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "discretePreOrderReport";
    private static final String TOLERANCE = "tolerance";

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
        return new DiscretePreOrderReport(likelihood, tolerance);
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
                AttributeRule.newDoubleRule(TOLERANCE, true)
        };
    }

    public static final XMLObjectParser PARSER = new DiscretePreOrderReportParser();
}
