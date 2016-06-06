package dr.inference.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;

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
    CompoundParameter traitsTemp;
    FullyConjugateMultivariateTraitLikelihood evaluator;


    public FactorRJMCMCOperator(CoercionMode mode, AdaptableSizeFastMatrixParameter factors, AdaptableSizeFastMatrixParameter loadings, AdaptableSizeFastMatrixParameter cutoffs, AdaptableSizeFastMatrixParameter loadingsSparcity, Distribution cutoffPrior, MomentDistributionModel loadingsPrior, FullyConjugateMultivariateTraitLikelihood factorsPrior) {
        super(mode);
        randomTree = new GaussianProcessFromTree(factorsPrior);
        this.factors = factors;
        this.loadings = loadings;
        this.cutoffs = cutoffs;
        this.loadingsSparcity = loadingsSparcity;
        this.cutoffPrior = cutoffPrior;
        this.loadingsPrior = loadingsPrior;
        Parameter[] paramListTemp = new Parameter.Default[1];
        paramListTemp[1] = new Parameter.Default(factors.getColumnDimension());
        this.traitsTemp = new CompoundParameter(null, paramListTemp);
        evaluator = factorsPrior.semiClone(traitsTemp);

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
        double from1 = 0;
        int currentSize = factors.getRowDimension();
        if(random > .5 || currentSize == 1){
            if(factors.getRowDimension() == 1) {
                from1 = Math.log(2);
            }
            return increment() + from1;
        }
        else{
            if(currentSize == 2){
                from1 = -Math.log(2);
            }
            return decrement() + from1;
        }



    }

    private double increment(){
        double tally = 0;
        int newSize = factors.getRowDimension() + 1;
        factors.setRowDimension(newSize);
        loadings.setColumnDimension(newSize);
        cutoffs.setColumnDimension(newSize);
        loadingsSparcity.setColumnDimension(newSize);


        //draws values from tree prior
        double[] temp = randomTree.nextRandomFast();
        for (int i = 0; i < temp.length; i++) {
            traitsTemp.setParameterValueQuietly(i, temp[i]);
        }
        traitsTemp.fireParameterChangedEvent();
        //computes likelihood of those values
        tally += evaluator.getLogLikelihood();

        //Draws new cutoff values
        for (int i = 0; i < cutoffs.getRowDimension(); i++) {
            double cutoffPriorDraw = cutoffPrior.quantile(MathUtils.nextDouble());
            tally += cutoffPrior.logPdf(cutoffPriorDraw);
        }

        //split apart rows
        int rowSplit = MathUtils.nextInt(newSize - 1);
        tally += Math.log((newSize - 1));
        for (int i = 0; i < loadingsSparcity.getRowDimension(); i++) {
            double splitDecider = MathUtils.nextDouble();
            if(splitDecider < .4){
                tally += Math.log(.4);
                loadingsSparcity.setParameterValueQuietly(i, newSize - 1, 0);
            }
            else if(splitDecider > .6){
                tally += Math.log(.4);
                loadingsSparcity.setParameterValueQuietly(i, rowSplit, 0);
                loadingsSparcity.setParameterValueQuietly(i, newSize - 1, 1);
            }
            else{
                tally += Math.log(.2);
                loadingsSparcity.setParameterValueQuietly(i, newSize - 1 , 1);
            }
        }


        for (int i = 0; i < loadings.getRowDimension(); i++) {
            NormalDistribution useful = getLoadingsDistribution(i, newSize - 1);
            double low = useful.cdf(-Math.sqrt(cutoffs.getParameterValue(i, newSize - 1)));
            double high = useful.cdf(Math.sqrt(cutoffs.getParameterValue(i, newSize - 1)));
            double proportion = low/(low + 1 - high);
            if(MathUtils.nextDouble() < proportion){
                double quantile=MathUtils.nextDouble() * low;
                loadings.setParameterValueQuietly(i, newSize - 1, useful.quantile(quantile));
                tally += useful.logPdf(loadings.getParameterValue(i, newSize - 1)) - Math.log(high - low);
            }
            else{
                double quantile=(1-high) * MathUtils.nextDouble() + high;
                loadings.setParameterValue(i, newSize - 1, useful.quantile(quantile));
                tally += useful.logPdf(loadings.getParameterValue(i, newSize - 1)) - Math.log(high - low);
            }

        }
        loadings.fireParameterChangedEvent();

        return tally;
    }

    private double decrement(){
        double tally = 0;
        int condensing1 = MathUtils.nextInt(factors.getRowDimension());
        int condensing2 = MathUtils.nextInt(factors.getRowDimension() - 1);
        tally -= Math.log(factors.getRowDimension() -1);

        if (condensing2 >= condensing1){
            condensing2 += 1;
        }
        
        
        for (int i = 0; i < loadings.getRowDimension(); i++) {
            NormalDistribution loadingsCondensed = getLoadingsDistribution(i, condensing2);
            double low = loadingsCondensed.cdf(-Math.sqrt(cutoffs.getParameterValue(i, condensing2)));
            double high = loadingsCondensed.cdf(Math.sqrt(cutoffs.getParameterValue(i, condensing2)));
            tally -= loadingsCondensed.logPdf(loadings.getParameterValue(i, condensing2)) - Math.log(high - low);
        }

        for (int i = 0; i < cutoffs.getRowDimension(); i++) {
            tally -= cutoffPrior.logPdf(cutoffs.getParameterValue(i, condensing2));
        }

        for (int i = 0; i < factors.getColumnDimension(); i++) {
            traitsTemp.setParameterValueQuietly(i, factors.getParameterValue(condensing2, i));
        }
        traitsTemp.fireParameterChangedEvent();
        tally -= evaluator.getLogLikelihood();


        for (int i = 0; i < loadingsSparcity.getColumnDimension(); i++) {
            if(loadingsSparcity.getParameterValue(i, condensing1) != loadingsSparcity.getParameterValue(i, condensing2))
            {
                tally -= Math.log(.4);
                loadingsSparcity.setParameterValueQuietly(i, condensing1, 1);
            }
            else if(loadingsSparcity.getParameterValue(i, condensing1) == 0){}
            else{
                tally -= Math.log(.2);
            }
        }

        for (int i = 0; i < factors.getRowDimension() ; i++) {
            factors.setParameterValueQuietly(condensing2, i, factors.getParameterValue(factors.getRowDimension()-1, i));
            loadings.setParameterValueQuietly(i, condensing2, loadings.getParameterValue(i, loadings.getColumnDimension()-1));
            loadingsSparcity.setParameterValueQuietly(i, condensing2, loadingsSparcity.getParameterValue(i, loadingsSparcity.getColumnDimension()-1));
            cutoffs.setParameterValueQuietly(i, condensing2, cutoffs.getParameterValue(i, cutoffs.getColumnDimension()-1));
        }


        return tally;
    }

    private NormalDistribution getLoadingsDistribution(int row, int column){
        //TODO fix LoadingsGibbsTruncatedOperator and apply it here
        double mean = 0;
        double variance = 1;
        return new NormalDistribution(mean, variance);
    }
}
