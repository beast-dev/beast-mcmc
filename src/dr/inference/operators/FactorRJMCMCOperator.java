package dr.inference.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;

/**
 * Created by max on 4/29/16.
 */
public class FactorRJMCMCOperator  extends AbstractCoercableOperator{
    GaussianProcessFromTree randomTree;
    FullyConjugateMultivariateTraitLikelihood factorsPrior;
    AdaptableSizeFastMatrixParameter factors;
    AdaptableSizeFastMatrixParameter loadings;
    AdaptableSizeFastMatrixParameter cutoffs;
    AdaptableSizeFastMatrixParameter loadingsSparcity;
    Distribution cutoffPrior;
    MomentDistributionModel loadingsPrior;


    public FactorRJMCMCOperator(CoercionMode mode, AdaptableSizeFastMatrixParameter factors, AdaptableSizeFastMatrixParameter loadings, AdaptableSizeFastMatrixParameter cutoffs, AdaptableSizeFastMatrixParameter loadingsSparcity, Distribution cutoffPrior, MomentDistributionModel loadingsPrior, FullyConjugateMultivariateTraitLikelihood factorsPrior) {
        super(mode);
        randomTree = new GaussianProcessFromTree(factorsPrior);
        this.factors = factors;
        this.loadings = loadings;
        this.cutoffs = cutoffs;
        this.loadingsSparcity = loadingsSparcity;
        this.cutoffPrior = cutoffPrior;
        this.loadingsPrior = loadingsPrior;

    }

    @Override
    public double getCoercableParameter() {
        return 0;
    }

    @Override
    public void setCoercableParameter(double value) {

    }

    @Override
    public double getRawParameter() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "FactorRJMCMCOperator";
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        double random = MathUtils.nextDouble();
        if(random > .5){
            if(factors.getRowDimension() == 1){}
        }
        Parameter.Default temp = new Parameter.Default(randomTree.nextRandomFast());


        return 0;
    }
}
