package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.treelikelihood.thorneytreelikelihood.*;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

public class ImportanceTreeSamplerParser extends AbstractXMLObjectParser {
    public static final String IMPORTANCE_TREE_SAMPLER = "importanceTreeSampler";
    public static final String SAMPLES = "samples";
    public static final String THREADS = "threads";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ThorneyTreeLikelihood treeLikelihood = (ThorneyTreeLikelihood) xo.getChild(ThorneyTreeLikelihood.class);

        GMRFMultilocusSkyrideLikelihood skygrid = (GMRFMultilocusSkyrideLikelihood) xo.getChild(GMRFMultilocusSkyrideLikelihood.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final int samples = xo.getAttribute(SAMPLES,1);
        final int threads = xo.getAttribute(THREADS,1);

        return new ImportanceTreeSampler(weight, samples, threads,treeLikelihood, skygrid);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "An operator that simulates a full tree from a constraints tree given the coalescent model and number of " +
                "mutations on each branch.";
    }

    @Override
    public Class getReturnType() {
        return AbstractTreeOperator.class;
    }

    @Override
    public String getParserName() {
        return IMPORTANCE_TREE_SAMPLER;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(SAMPLES, true),
            AttributeRule.newIntegerRule(THREADS, true),
            new ElementRule(ThorneyTreeLikelihood.class),
            new ElementRule(GMRFMultilocusSkyrideLikelihood.class),
    };
}
