package dr.inference.operators;

import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.distribution.RowDimensionPoissonPrior;
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
    RowDimensionPoissonPrior rowPrior;

    public static final boolean DEBUG = false;

    public FactorRJMCMCOperator(double weight, double sizeParam, int chainLength, AdaptableSizeFastMatrixParameter factors, AdaptableSizeFastMatrixParameter loadings, AdaptableSizeFastMatrixParameter cutoffs, AdaptableSizeFastMatrixParameter loadingsSparsity, LatentFactorModel lfm, DeterminentalPointProcessPrior sparsityPrior, LoadingsGibbsTruncatedOperator loadingsOperator, FactorTreeGibbsOperator factorOperator, BitFlipOperator sparsityOperator, RowDimensionPoissonPrior rowPrior) {
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
        this.rowPrior = rowPrior;
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

        if (DEBUG) {
            System.out.println("sparsity prior lnL = " + sparsityPrior.getLogLikelihood());
            System.out.println("latentFactorModel lnL = " + lfm.getLogLikelihood());
            System.out.print("loadingsSparsity: ");
            for (int i = 0; i <loadingsSparsity.getDimension() ; i++) {
                System.out.print(loadingsSparsity.getParameterValue(i) + " ");
            }
        }

        double random = MathUtils.nextDouble();
        double from1 = 0;
        double initialLikelihood = lfm.getLogLikelihood() * (1 - sizeParam) + rowPrior.getLogLikelihood();
        boolean increment;

        if (DEBUG) {
            System.out.print("\nloadings: ");
            for (int i = 0; i < loadings.getDimension() ; i++) {
                System.out.print(loadings.getParameterValue(i) + " ");
            }
        }

        storeDimensions();
        storeValues();

        if (DEBUG) {
            System.out.print("\nstoredLoadings: ");
            for (int i = 0; i < storedLoadings.getDimension(); i++) {
                System.out.print(storedLoadings.getParameterValue(i) + " ");
            }
            System.out.println("\nstoredLoadingsSparsity: ");
            for (int i = 0; i < storedLoadingsSparsity.getDimension(); i++) {
                System.out.print(storedLoadingsSparsity.getParameterValue(i) + " ");
            }
            System.out.println();
        }

        int currentSize = factors.getRowDimension();
        if((random > .5 || currentSize == 1) && currentSize != factors.getMaxRowDimension()){
            if(factors.getRowDimension() == 1) {
                from1 = Math.log(2);
            }
            if(factors.getRowDimension() == factors.getRowDimension() - 1){
                from1 = -Math.log(2);
            }
            factors.setRowDimension(factors.getRowDimension()+1);
            loadings.setColumnDimension(loadings.getColumnDimension()+1);
            loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()+1);
            cutoffs.setColumnDimension(cutoffs.getColumnDimension()+1);
            if (DEBUG) {
                System.out.println("up");
            }
            increment = true;
        }
        else{
            if(factors.getRowDimension() == factors.getMaxRowDimension()) {
                from1 = Math.log(2);
            }
            if(currentSize == 2 || currentSize == factors.getDimension()){
                from1 = -Math.log(2);
            }
            factors.setRowDimension(factors.getRowDimension()-1);
            loadings.setColumnDimension(loadings.getColumnDimension()-1);
            loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()-1);
            cutoffs.setColumnDimension(cutoffs.getColumnDimension()-1);
            if (DEBUG) {
                System.out.println("down");
            }
            increment = false;
        }

        //hack to let me store model state later in the code
        lfm.acceptModelState();
        sparsityPrior.acceptModelState();

        iterate();

        if (DEBUG) {
            System.out.print("storedFactors: ");
            for (int i = 0; i <storedFactors.getDimension(); i++) {
                System.out.print(storedFactors.getParameterValue(i) + " ");
            }
            System.out.println();
        }

        double finalLikelihood = lfm.getLogLikelihood() * (1 - sizeParam) + rowPrior.getLogLikelihood();

        random = MathUtils.nextDouble();
        double test = from1 + finalLikelihood - initialLikelihood;
        test = Math.min(Math.exp(test), 1);
        if(random < test){
            if (DEBUG) {
                System.out.println("accepted!\n" + test);
            }
            lfm.acceptModelState();
            lfm.makeDirty();
            sparsityPrior.acceptModelState();
        }
        else{
            restoreDimensions();
            restoreValues();
            if (DEBUG) {
                System.out.print("loadingsSparsity: ");
                for (int i = 0; i <loadingsSparsity.getDimension(); i++) {
                    System.out.print(loadingsSparsity.getParameterValue(i) + " ");
                }
                System.out.println();
            }
            sparsityPrior.makeDirty();
            lfm.makeDirty();
            if (DEBUG) {
                System.out.println("sparsity prior lnL = " + sparsityPrior.getLogLikelihood());
                System.out.println("latentFactorModel lnL = " + lfm.getLogLikelihood());
                System.out.print("loadingsSparsity: ");
                for (int i = 0; i <loadingsSparsity.getDimension(); i++) {
                    System.out.print(loadingsSparsity.getParameterValue(i) + " ");
                }
                System.out.println();
            }
            factors.storeParameterValues();
            loadings.storeParameterValues();
            loadingsSparsity.storeParameterValues();
            cutoffs.storeParameterValues();
            if (DEBUG) {
                System.out.println("rejected!\n" + test);
            }

        }
        return 0;
    }

    private void iterate() throws OperatorFailedException {
        factorOperator.setPathParameter(sizeParam);
        loadingsOperator.setPathParameter(sizeParam);
//        if(separator == null){
            separator = new double[2];
            double foWeight = factors.getColumnDimension() * chainLength;
            double loWeight = loadings.getColumnDimension() * chainLength;
            double sparoWeight = (loadings.getRowDimension() * loadings.getColumnDimension() + factors.getColumnDimension() * factors.getRowDimension()) * chainLength;
            double total = foWeight + loWeight + sparoWeight;
            separator[0] = foWeight / total;
            separator[1] = (foWeight + loWeight) / total;
//        }
        for (int i = 0; i < total; i++) {
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
                mhRatio += sparsityOperator.doOperation();
                mhRatio += lfm.getLogLikelihood() * sizeParam + sparsityPrior.getLogLikelihood();
                mhRatio = Math.min(1, Math.exp(mhRatio));
                if(MathUtils.nextDouble() > mhRatio || Double.isNaN(sparsityPrior.getLogLikelihood())){
                    lfm.restoreModelState();
                    sparsityPrior.restoreModelState();
//                    lfm.makeDirty();
//                    sparsityPrior.makeDirty();
                }
                else{
                    lfm.acceptModelState();
                    sparsityPrior.acceptModelState();
//                    lfm.makeDirty();
//                    sparsityPrior.makeDirty();
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
        storedLoadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension());
        storedCutoffs.setRowDimension(cutoffs.getRowDimension());
        storedCutoffs.setColumnDimension(cutoffs.getColumnDimension());}

    private void restoreDimensions(){
        factors.setRowDimension(storedFactors.getRowDimension());
        factors.setColumnDimension(storedFactors.getColumnDimension());
        loadings.setRowDimension(storedLoadings.getRowDimension());
        loadings.setColumnDimension(storedLoadings.getColumnDimension());
        loadingsSparsity.setRowDimension(storedLoadingsSparsity.getRowDimension());
        loadingsSparsity.setColumnDimension(storedLoadingsSparsity.getColumnDimension());
        cutoffs.setRowDimension(storedCutoffs.getRowDimension());
        cutoffs.setColumnDimension(storedCutoffs.getColumnDimension());}


    private void storeValues(){
        for (int i = 0; i < factors.getDimension(); i++) {
            storedFactors.setParameterValue(i, factors.getParameterValue(i));}
        for (int i = 0; i < loadings.getDimension(); i++) {
            storedLoadings.setParameterValue(i, loadings.getParameterValue(i));
            storedLoadingsSparsity.setParameterValue(i, loadingsSparsity.getParameterValue(i));
            storedCutoffs.setParameterValue(i, cutoffs.getParameterValue(i));
        }
    }

    private void restoreValues(){
        for (int i = 0; i < factors.getDimension(); i++) {
            factors.setParameterValue(i, storedFactors.getParameterValue(i));}
        for (int i = 0; i < loadings.getDimension(); i++) {
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
