package dr.inferencexml.operators;

import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.BitFlipOperator;
import dr.inference.operators.FactorRJMCMCOperator;
import dr.inference.operators.FactorTreeGibbsOperator;
import dr.inference.operators.LoadingsGibbsTruncatedOperator;
import dr.xml.*;

/**
 * Created by maxryandolinskytolkoff on 6/13/16.
 */
public class FactorRJMCMCOperatorParser extends AbstractXMLObjectParser{
    public static final String FACTOR_RJMCMC_OPERATOR = "factorRJMCMCOperator";
    public static final String WEIGHT = "weight";
    public static final String LOADINGS_SPARSITY = "loadingsSparsity";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String FACTORS = "factors";
    public static final String LOADINGS = "loadings";
    public static final String CUTOFFS = "cutoffs";
    public static final String SIZE_PARAMETER = "sizeParameter";



    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = xo.getDoubleAttribute(WEIGHT);
        int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);
        AdaptableSizeFastMatrixParameter factors, loadings, cutoffs, loadingsSparcity;
        factors = (AdaptableSizeFastMatrixParameter) xo.getChild(FACTORS).getChild(AdaptableSizeFastMatrixParameter.class);
        loadings = (AdaptableSizeFastMatrixParameter) xo.getChild(LOADINGS).getChild(AdaptableSizeFastMatrixParameter.class);
        cutoffs = (AdaptableSizeFastMatrixParameter) xo.getChild(CUTOFFS).getChild(AdaptableSizeFastMatrixParameter.class);
        loadingsSparcity = (AdaptableSizeFastMatrixParameter) xo.getChild(LOADINGS_SPARSITY).getChild(AdaptableSizeFastMatrixParameter.class);
        DeterminentalPointProcessPrior DPP = (DeterminentalPointProcessPrior) xo.getChild(DeterminentalPointProcessPrior.class);
        LatentFactorModel LFM = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        BitFlipOperator sparsityOperator = (BitFlipOperator) xo.getChild(BitFlipOperator.class);
        LoadingsGibbsTruncatedOperator loadingsOperator = (LoadingsGibbsTruncatedOperator) xo.getChild(LoadingsGibbsTruncatedOperator.class);
        FactorTreeGibbsOperator factorOperator = (FactorTreeGibbsOperator) xo.getChild(FactorTreeGibbsOperator.class);
        double sizeParameter = xo.getDoubleAttribute(SIZE_PARAMETER);



        return new FactorRJMCMCOperator(weight, sizeParameter, chainLength, factors, loadings, cutoffs, loadingsSparcity, LFM, DPP, loadingsOperator, factorOperator, sparsityOperator);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LoadingsGibbsTruncatedOperator.class),
            new ElementRule(FactorTreeGibbsOperator.class),
            new ElementRule(BitFlipOperator.class),
            new ElementRule(LatentFactorModel.class),
            new ElementRule(DeterminentalPointProcessPrior.class),
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(CUTOFFS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(LOADINGS_SPARSITY, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newIntegerRule(CHAIN_LENGTH),
            AttributeRule.newDoubleRule(SIZE_PARAMETER),
    };

    @Override
    public String getParserDescription() {
        return "RJMCMC to determine the number of factors in a factor analysis model";
    }

    @Override
    public Class getReturnType() {
        return FactorRJMCMCOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_RJMCMC_OPERATOR;
    }
}
