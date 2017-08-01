package dr.inferencexml.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.operators.LFMSplitMergeOperator;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class LFMSplitMergeOperatorParser extends AbstractXMLObjectParser {
    public final static String LFM_SPLIT_MERGE_OPERATOR = "LFMSplitMergeOperator";
    public final static String FACTORS = "factors";
    public final static String LOADINGS = "loadings";
    public final static String CUTOFFS = "cutoffs";
    public final static String SPARSE_LOADINGS = "sparseLoadings";
    public static final String WEIGHT = "weight";
    public static final String SPLIT_VARIANCE = "splitVariance";

    @Override
    public String getParserName() {
        return LFM_SPLIT_MERGE_OPERATOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AdaptableSizeFastMatrixParameter loadings = (AdaptableSizeFastMatrixParameter) xo.getChild(LOADINGS).getChild(AdaptableSizeFastMatrixParameter.class);
        AdaptableSizeFastMatrixParameter factors = (AdaptableSizeFastMatrixParameter) xo.getChild(FACTORS).getChild(AdaptableSizeFastMatrixParameter.class);
        AdaptableSizeFastMatrixParameter cutoffs = (AdaptableSizeFastMatrixParameter) xo.getChild(CUTOFFS).getChild(AdaptableSizeFastMatrixParameter.class);
        AdaptableSizeFastMatrixParameter sparseLoadings = (AdaptableSizeFastMatrixParameter) xo.getChild(SPARSE_LOADINGS).getChild(AdaptableSizeFastMatrixParameter.class);
        double weight = xo.getDoubleAttribute(WEIGHT);
        double splitVariance = xo.getDoubleAttribute(SPLIT_VARIANCE);
        FullyConjugateMultivariateTraitLikelihood tree = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

        return new LFMSplitMergeOperator(weight, splitVariance, factors, loadings, sparseLoadings, cutoffs, tree);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }


    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(SPARSE_LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(CUTOFFS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(SPLIT_VARIANCE),
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class)
    };







    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return LFMSplitMergeOperator.class;
    }
}
