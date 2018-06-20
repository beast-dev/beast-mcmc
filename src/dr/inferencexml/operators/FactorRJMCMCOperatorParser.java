package dr.inferencexml.operators;

import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.distribution.RowDimensionPoissonPrior;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.Likelihood;
import dr.inference.operators.*;
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
    public static final String LOADINGS_PRIOR = "loadingsPrior";
    public static final String FACTOR_OPERATOR = "factorOperator";
    public static final String LOADINGS_OPERATOR = "loadingsOperator";
    public static final String ROW_PRIOR = "rowPrior";
    public static final String SPARSITY_PRIOR = "sparsityPrior";
    public static final String NEGATION_OPERATOR = "negationOperator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        //attributes
        double weight = xo.getDoubleAttribute(WEIGHT);
        int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);
        double sizeParameter = xo.getDoubleAttribute(SIZE_PARAMETER);

        //declaration
        AdaptableSizeFastMatrixParameter factors, loadings, cutoffs, loadingsSparcity;


        //parameters
        if(xo.hasChildNamed(FACTORS))
            factors = (AdaptableSizeFastMatrixParameter) xo.getChild(FACTORS).getChild(AdaptableSizeFastMatrixParameter.class);
        else
            factors = null;


        loadings = (AdaptableSizeFastMatrixParameter) xo.getChild(LOADINGS).getChild(AdaptableSizeFastMatrixParameter.class);

        if(xo.hasChildNamed(CUTOFFS)){
            cutoffs = (AdaptableSizeFastMatrixParameter) xo.getChild(CUTOFFS).getChild(AdaptableSizeFastMatrixParameter.class);}
        else
            cutoffs = null;
        if(xo.hasChildNamed(LOADINGS_SPARSITY))
            loadingsSparcity = (AdaptableSizeFastMatrixParameter) xo.getChild(LOADINGS_SPARSITY).getChild(AdaptableSizeFastMatrixParameter.class);
        else
            loadingsSparcity = null;

        //models
        DeterminentalPointProcessPrior DPP = null;
        if(xo.getChild(SPARSITY_PRIOR) != null)
             DPP = (DeterminentalPointProcessPrior) xo.getChild(SPARSITY_PRIOR).getChild(DeterminentalPointProcessPrior.class);
        SimpleMCMCOperator NOp= null;
        if(xo.getChild(NEGATION_OPERATOR) != null){
            NOp = (SimpleMCMCOperator) xo.getChild(NEGATION_OPERATOR).getChild(SimpleMCMCOperator.class);
        }
        AbstractModelLikelihood LFM = (AbstractModelLikelihood) xo.getChild(AbstractModelLikelihood.class);
        RowDimensionPoissonPrior rowPrior = (RowDimensionPoissonPrior) xo.getChild(ROW_PRIOR).getChild(RowDimensionPoissonPrior.class);
        Likelihood loadingsPrior = null;
        if(xo.hasChildNamed(LOADINGS_PRIOR)){
            loadingsPrior = (Likelihood) xo.getChild(LOADINGS_PRIOR).getChild(Likelihood.class);
        }

        //operators
        SimpleMCMCOperator sparsityOperator = null;
        if (xo.getChild(BitFlipOperator.class) != null)
            sparsityOperator = (BitFlipOperator) xo.getChild(BitFlipOperator.class);
        if(xo.getChild(LoadingsSparsityOperator.class) != null){
            sparsityOperator = (LoadingsSparsityOperator) xo.getChild(LoadingsSparsityOperator.class);
        }
        SimpleMCMCOperator loadingsOperator = (SimpleMCMCOperator) xo.getChild(LOADINGS_OPERATOR).getChild(SimpleMCMCOperator.class);
        SimpleMCMCOperator factorOperator = null;
        if(xo.getChild(FACTOR_OPERATOR) != null)
            factorOperator = (SimpleMCMCOperator) xo.getChild(FACTOR_OPERATOR).getChild(FactorTreeGibbsOperator.class);

        LatentFactorModelPrecisionGibbsOperator precisionGibbsOperator = null;
        if(xo.getChild(LatentFactorModelPrecisionGibbsOperator.class) != null){
            System.out.println("here");
            precisionGibbsOperator = (LatentFactorModelPrecisionGibbsOperator) xo.getChild(LatentFactorModelPrecisionGibbsOperator.class);
        }




        return new FactorRJMCMCOperator(weight, sizeParameter, chainLength, factors,
                loadings, cutoffs, loadingsSparcity, LFM, DPP,
                loadingsPrior, loadingsOperator, factorOperator,
                sparsityOperator, NOp, rowPrior, precisionGibbsOperator);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOADINGS_OPERATOR, new XMLSyntaxRule[]{
                    new ElementRule(SimpleMCMCOperator.class)}),
            new ElementRule(FACTOR_OPERATOR, new XMLSyntaxRule[]{
                    new ElementRule(SimpleMCMCOperator.class)}, true),
            new OrRule(new ElementRule(BitFlipOperator.class, true), new ElementRule(LoadingsSparsityOperator.class, true)),
            new ElementRule(AbstractModelLikelihood.class),
            new ElementRule(LatentFactorModelPrecisionGibbsOperator.class, true),
            new ElementRule(SPARSITY_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(DeterminentalPointProcessPrior.class)}, true),
            new ElementRule(FACTORS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}, true),
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}),
            new ElementRule(CUTOFFS, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}, true),
            new ElementRule(LOADINGS_SPARSITY, new XMLSyntaxRule[]{
                    new ElementRule(AdaptableSizeFastMatrixParameter.class)}, true),
            new ElementRule(NEGATION_OPERATOR, new XMLSyntaxRule[]{
                    new ElementRule(SimpleMCMCOperator.class)}, true),
            new ElementRule(ROW_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(RowDimensionPoissonPrior.class)}),
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
