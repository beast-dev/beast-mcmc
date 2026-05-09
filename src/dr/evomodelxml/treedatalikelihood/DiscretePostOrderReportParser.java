
package dr.evomodelxml.treedatalikelihood;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.DiscretePostOrderReport;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Reportable XML helper for dumping all cached discrete post-order partials.
 */
public class DiscretePostOrderReportParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "discretePostOrderReport";

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
        return new DiscretePostOrderReport(likelihood);
    }

    @Override
    public String getParserDescription() {
        return "Reports branch-start and branch-end post-order partials for a discrete TreeDataLikelihood.";
    }

    @Override
    public Class getReturnType() {
        return DiscretePostOrderReport.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    public static final XMLObjectParser PARSER = new DiscretePostOrderReportParser();
}
