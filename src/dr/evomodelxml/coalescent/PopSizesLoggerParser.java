package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.GMRFSkygridLikelihood;
import dr.evomodel.coalescent.MultilocusNonparametricCoalescentLikelihood;
import dr.evomodel.coalescent.PopSizesLogger;
import dr.inference.model.Parameter;
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
        if (xo.getChild(MultilocusNonparametricCoalescentLikelihood.class) != null) {
            MultilocusNonparametricCoalescentLikelihood likelihood = (MultilocusNonparametricCoalescentLikelihood) xo.getChild(MultilocusNonparametricCoalescentLikelihood.class);
            return new PopSizesLogger(likelihood, transform, order);
        } else if (xo.getChild(GMRFSkygridLikelihood.class) != null) {
            GMRFSkygridLikelihood likelihood = (GMRFSkygridLikelihood) xo.getChild(GMRFSkygridLikelihood.class);
            return new PopSizesLogger(likelihood, transform, order);
        } else if (xo.getChild(GMRFMultilocusSkyrideLikelihood.class) != null) {
            GMRFMultilocusSkyrideLikelihood likelihood = (GMRFMultilocusSkyrideLikelihood) xo.getChild(GMRFMultilocusSkyrideLikelihood.class);
            return new PopSizesLogger(likelihood, transform, order);
        } else {
            throw new XMLParseException("Likelihood type not yet implemented.");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(ORDER, true),
            new XORRule(
                    new ElementRule(MultilocusNonparametricCoalescentLikelihood.class),
                    new XORRule(
                            new ElementRule(GMRFSkygridLikelihood.class),
                            new ElementRule(GMRFMultilocusSkyrideLikelihood.class))
            ),

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
