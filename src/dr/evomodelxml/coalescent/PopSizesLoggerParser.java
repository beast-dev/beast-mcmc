package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.MultilocusNonparametricCoalescentLikelihood;
import dr.evomodel.coalescent.PopSizesLogger;
import dr.util.Transform;
import dr.xml.*;

public class PopSizesLoggerParser extends AbstractXMLObjectParser {

    private static final String NAME = "popSizesLogger";
    private static final String ORDER = "order";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String order = xo.getAttribute(ORDER, "backwards");
        Transform.ParsedTransform pt = (Transform.ParsedTransform) xo.getChild(Transform.ParsedTransform.class);
        Transform transform = null;
        if (pt != null) {
            transform = pt.transform;
        }
        MultilocusNonparametricCoalescentLikelihood likelihood = (MultilocusNonparametricCoalescentLikelihood) xo.getChild(MultilocusNonparametricCoalescentLikelihood.class);
        return new PopSizesLogger(likelihood, transform, order);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(ORDER, true),
            new ElementRule(MultilocusNonparametricCoalescentLikelihood.class),
            new ElementRule(Transform.ParsedTransform.class, true)
    };

    @Override
    public String getParserDescription() {
        return "Logger to report population sizes of a coalescent model";
    }

    @Override
    public Class getReturnType() {
        return PopSizesLogger.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
