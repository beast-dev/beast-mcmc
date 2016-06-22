package dr.inference.operators;

import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

/**
 * Created by max on 4/29/16.
 */
public class FactorRJMCMCOperator  extends SimpleMCMCOperator implements GibbsOperator{
    GaussianProcessFromTree randomTree;
    AdaptableSizeFastMatrixParameter factors;
    AdaptableSizeFastMatrixParameter loadings;
    AdaptableSizeFastMatrixParameter cutoffs;
    AdaptableSizeFastMatrixParameter loadingsSparsity;
    LatentFactorModel lfm;
    DeterminentalPointProcessPrior sparsityPrior;
    int chainLength;
    CompoundParameter traitsTemp;
    double sizeParam;
    private double[] separator;
    LoadingsGibbsTruncatedOperator loadingsOperator;
    FactorTreeGibbsOperator factorOperator;
    BitFlipOperator sparsityOperator;
    AdaptableSizeFastMatrixParameter storedFactors;
    AdaptableSizeFastMatrixParameter storedLoadings;
    AdaptableSizeFastMatrixParameter storedCutoffs;
    AdaptableSizeFastMatrixParameter storedLoadingsSparsity;


    public FactorRJMCMCOperator(double weight, double sizeParam, int chainLength, AdaptableSizeFastMatrixParameter factors, AdaptableSizeFastMatrixParameter loadings, AdaptableSizeFastMatrixParameter cutoffs, AdaptableSizeFastMatrixParameter loadingsSparsity,  LatentFactorModel lfm, DeterminentalPointProcessPrior sparsityPrior, LoadingsGibbsTruncatedOperator loadingsOperator, FactorTreeGibbsOperator factorOperator, BitFlipOperator sparsityOperator) {
        setWeight(weight);
        this.factors = factors;
        this.loadings = loadings;
        this.cutoffs = cutoffs;
        this.loadingsSparsity = loadingsSparsity;
        this.sparsityPrior = sparsityPrior;
        this.lfm = lfm;
        this.sizeParam = sizeParam;
        this.chainLength = chainLength;
//        Parameter[] paramListTemp = new Parameter.Default[1];
//        paramListTemp[1] = new Parameter.Default(factors.getColumnDimension());
//        this.traitsTemp = new CompoundParameter(null, paramListTemp);
        this.storedFactors = new AdaptableSizeFastMatrixParameter(factors.getId()+".stored", 1, 1, factors.getMaxRowDimension(), factors.getMaxColumnDimension(), 1);
        this.storedLoadings = new AdaptableSizeFastMatrixParameter(loadings.getId()+".stored", 1, 1, loadings.getMaxRowDimension(), loadings.getMaxColumnDimension(), 1);
        this.storedCutoffs = new AdaptableSizeFastMatrixParameter(cutoffs.getId()+".stored", 1, 1, cutoffs.getMaxRowDimension(), cutoffs.getMaxColumnDimension(), 1);
        this.storedLoadingsSparsity = new AdaptableSizeFastMatrixParameter(loadingsSparsity.getId()+".stored", 1, 1, loadingsSparsity.getMaxRowDimension(), loadingsSparsity.getMaxColumnDimension(), 1);
        this.loadingsOperator = loadingsOperator;
        this.factorOperator = factorOperator;
        this.sparsityOperator = sparsityOperator;
        storeDimensions();
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
        double initialLikelihood = lfm.getLogLikelihood() * (1 - sizeParam);
        boolean increment;



        int currentSize = factors.getRowDimension();
        if(random > .5 || currentSize == 1){
            if(factors.getRowDimension() == 1) {
                from1 = Math.log(2);
            }
            storeDimensions();
            storeValues();
            factors.setRowDimension(factors.getRowDimension()+1);
            loadings.setColumnDimension(loadings.getColumnDimension()+1);
            loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()+1);
            cutoffs.setColumnDimension(cutoffs.getColumnDimension()+1);
            System.out.println("up");
            increment = true;
        }
        else{
            if(currentSize == 2){
                from1 = -Math.log(2);
            }
            storeValues();
            storeDimensions();
            factors.setRowDimension(factors.getRowDimension()-1);
            loadings.setColumnDimension(loadings.getColumnDimension()-1);
            loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()-1);
            cutoffs.setColumnDimension(cutoffs.getColumnDimension()-1);
            System.out.println("down");
            increment = false;
        }

        iterate();
        double finalLikelihood = lfm.getLogLikelihood() *(1 - sizeParam);


        random = MathUtils.nextDouble();
        double test = from1 + finalLikelihood - initialLikelihood;
        test = Math.min(Math.exp(test), 1);
        if(random < test){
            System.out.println("yup!");
            System.out.println(test);
        }
        else{
            if(increment){
                restoreValues();
                restoreDimensions();
            }
            else{
                restoreDimensions();
                restoreValues();
            }

            sparsityPrior.makeDirty();
            lfm.makeDirty();
            factors.storeParameterValues();
            loadings.storeParameterValues();
            loadingsSparsity.storeParameterValues();
            cutoffs.storeParameterValues();
            System.out.println("nope :(");
            System.out.println(test);

        }
        return 0;
    }

    private void iterate() throws OperatorFailedException {
        factorOperator.setPathParameter(sizeParam);
        loadingsOperator.setPathParameter(sizeParam);
        if(separator == null){
            separator = new double[2];
            double foWeight = factorOperator.getWeight();
            double loWeight = loadingsOperator.getWeight();
            double sparoWeight = sparsityOperator.getWeight();
            double total = foWeight + loWeight + sparoWeight;
            separator[0] = foWeight/total;
            separator[1] = foWeight + loWeight;
        }
        for (int i = 0; i < chainLength; i++) {
            double rand = MathUtils.nextDouble();
            if(rand < separator[0]){
                factorOperator.doOperation();
            }
            else if (rand < separator[1]){
                loadingsOperator.doOperation();
            }
            else {
                lfm.storeModelState();
                sparsityPrior.storeModelState();
                double mhRatio = - lfm.getLogLikelihood() * sizeParam - sparsityPrior.getLogLikelihood();
                mhRatio -= sparsityOperator.doOperation();
                mhRatio += lfm.getLogLikelihood() * sizeParam -sparsityPrior.getLogLikelihood();
                mhRatio = Math.min(1, Math.exp(mhRatio));
                if(MathUtils.nextDouble() > mhRatio){
                    lfm.restoreModelState();
                    sparsityPrior.restoreModelState();
                }
            }

        }
        factorOperator.setPathParameter(1);
        loadingsOperator.setPathParameter(1);
    }

    private void storeDimensions(){storedFactors.setRowDimension(factors.getRowDimension());
        storedFactors.setColumnDimension(factors.getColumnDimension());
        storedLoadings.setRowDimension(loadings.getRowDimension());
        storedLoadings.setColumnDimension(loadings.getColumnDimension());
        storedLoadingsSparsity.setRowDimension(loadingsSparsity.getRowDimension());
        storedLoadingsSparsity.setColumnDimension(loadings.getColumnDimension());
        storedCutoffs.setRowDimension(cutoffs.getRowDimension());
        storedCutoffs.setColumnDimension(cutoffs.getColumnDimension());}

    private void restoreDimensions(){factors.setRowDimension(storedFactors.getRowDimension());
        factors.setColumnDimension(storedFactors.getColumnDimension());
        loadings.setRowDimension(storedLoadings.getRowDimension());
        loadings.setColumnDimension(storedLoadings.getColumnDimension());
        loadingsSparsity.setRowDimension(storedLoadingsSparsity.getRowDimension());
        loadingsSparsity.setColumnDimension(storedLoadingsSparsity.getColumnDimension());
        cutoffs.setRowDimension(storedCutoffs.getRowDimension());
        cutoffs.setColumnDimension(storedCutoffs.getColumnDimension());}

    private void storeValues(){
        for (int i = 0; i < factors.getDimension(); i++) {
            storedFactors.setParameterValue(i, factors.getParameterValue(i));
            storedLoadings.setParameterValue(i, loadings.getParameterValue(i));
            storedLoadingsSparsity.setParameterValue(i, loadingsSparsity.getParameterValue(i));
            storedCutoffs.setParameterValue(i, cutoffs.getParameterValue(i));
        }
    }

    private void restoreValues(){
        for (int i = 0; i < factors.getDimension(); i++) {
            factors.setParameterValue(i, storedFactors.getParameterValue(i));
            loadings.setParameterValue(i, storedLoadings.getParameterValue(i));
            loadingsSparsity.setParameterValue(i, storedLoadingsSparsity.getParameterValue(i));
            cutoffs.setParameterValue(i, storedCutoffs.getParameterValue(i));
        }
    }

    @Override
    public int getStepCount() {
        return 0;
    }
}
