package dr.inference.operators;

import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.distribution.DeterminentalPointProcessPrior;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.*;
import dr.math.MathUtils;

/**
 * Created by max on 4/29/16.
 */
public class FactorRJMCMCOperator  extends SimpleMCMCOperator implements GibbsOperator{

    GaussianProcessFromTree randomTree;
    AdaptableSizeFastMatrixParameter factors;
    AdaptableSizeFastMatrixParameter loadings;
    AdaptableSizeFastMatrixParameter cutoffs;
    AdaptableSizeFastMatrixParameter loadingsSparsity;
    AbstractModelLikelihood lfm;
    DeterminentalPointProcessPrior sparsityPrior;
    Likelihood loadingsPrior;
    int chainLength;
    CompoundParameter traitsTemp;
    double sizeParam;
    private double[] separator;
    SimpleMCMCOperator loadingsOperator;
    SimpleMCMCOperator factorOperator;
    SimpleMCMCOperator sparsityOperator;
    SimpleMCMCOperator NOp;
    AdaptableSizeFastMatrixParameter storedFactors;
    AdaptableSizeFastMatrixParameter storedLoadings;
    AdaptableSizeFastMatrixParameter storedCutoffs;
    AdaptableSizeFastMatrixParameter storedLoadingsSparsity;
    MatrixSizePrior rowPrior;
    LatentFactorModelPrecisionGibbsOperator precisionGibbsOperator;

    private final int BASE_SIZE = 1000;
    private final double MIN_WEIGHT = .01;
    private int callCount = 0;
    double callWeighting = 1;
    public static final boolean DEBUG = false;


    public FactorRJMCMCOperator(double weight, double sizeParam, int chainLength, AdaptableSizeFastMatrixParameter factors,
                                AdaptableSizeFastMatrixParameter loadings, AdaptableSizeFastMatrixParameter cutoffs,
                                AdaptableSizeFastMatrixParameter loadingsSparsity, AbstractModelLikelihood lfm,
                                DeterminentalPointProcessPrior sparsityPrior, Likelihood loadingsPrior,
                                SimpleMCMCOperator loadingsOperator, SimpleMCMCOperator factorOperator,
                                SimpleMCMCOperator sparsityOperator, SimpleMCMCOperator NOp, MatrixSizePrior rowPrior,
                                LatentFactorModelPrecisionGibbsOperator precisionGibbsOperator) {
        setWeight(weight);
        this.factors = factors;
        this.loadings = loadings;
        this.cutoffs = cutoffs;
        this.loadingsSparsity = loadingsSparsity;
        this.sparsityPrior = sparsityPrior;
        this.lfm = lfm;
        this.sizeParam = sizeParam;
        this.chainLength = chainLength;
        this.NOp = NOp;
//        Parameter[] paramListTemp = new Parameter.Default[1];
//        paramListTemp[1] = new Parameter.Default(factors.getColumnDimension());
//        this.traitsTemp = new CompoundParameter(null, paramListTemp);
        if(factors != null)
            this.storedFactors = new AdaptableSizeFastMatrixParameter(factors.getId()+".stored", 1, 1, factors.getMaxRowDimension(), factors.getMaxColumnDimension(), 1, false);
        this.storedLoadings = new AdaptableSizeFastMatrixParameter(loadings.getId()+".stored", 1, 1, loadings.getMaxRowDimension(), loadings.getMaxColumnDimension(), 1, false);
        if(cutoffs != null)
            this.storedCutoffs = new AdaptableSizeFastMatrixParameter(cutoffs.getId()+".stored", 1, 1, cutoffs.getMaxRowDimension(), cutoffs.getMaxColumnDimension(), 1, false);
        if(loadingsSparsity != null)
            this.storedLoadingsSparsity = new AdaptableSizeFastMatrixParameter(loadingsSparsity.getId()+".stored", 1, 1, loadingsSparsity.getMaxRowDimension(), loadingsSparsity.getMaxColumnDimension(), 1, false);
        this.loadingsOperator = loadingsOperator;
        this.factorOperator = factorOperator;
        this.sparsityOperator = sparsityOperator;
        this.rowPrior = rowPrior;
        this.loadingsPrior = loadingsPrior;
        this.precisionGibbsOperator = precisionGibbsOperator;
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
    public double doOperation() {//throws OperatorFailedException {
        boolean off = true;
        if(callCount < BASE_SIZE || off == true){
            performOperation();
        }
        else{
            callWeighting *= .99;
            if(callWeighting < MIN_WEIGHT){
                callWeighting = MIN_WEIGHT;
            }
            if(callWeighting > MathUtils.nextDouble()){
                performOperation();
            }
        }

        callCount++;
        return 0;
    }

    private void performOperation(){
        //        System.out.println(sparsityPrior.getLogLikelihood());
//        System.out.println(lfm.getLogLikelihood());
        String outpu="";
//        for (int i = 0; i <loadingsSparsity.getDimension() ; i++) {
//            outpu += loadingsSparsity.getParameterValue(i);
//            outpu += " ";
//        }
//        System.out.println(outpu);

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


//        System.out.println("Before");
//        System.out.println(lfm.getLogLikelihood() * (1 - sizeParam));
//        System.out.println(rowPrior.getSizeLogLikelihood());


        double initialLikelihood = lfm.getLogLikelihood() * (1 - sizeParam) + rowPrior.getSizeLogLikelihood();
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

        int currentSize = loadings.getColumnDimension();
        if((random > .5 || currentSize == 1) && currentSize != loadings.getMaxColumnDimension()){
            if(loadings.getColumnDimension() == 1) {
                from1 = -Math.log(2);
            }
            if(loadings.getColumnDimension() == loadings.getMaxColumnDimension() - 1){
                from1 = Math.log(2);
            }
            if(factors != null)
                factors.setRowDimension(factors.getRowDimension()+1);
            loadings.setColumnDimension(loadings.getColumnDimension()+1);
            if(loadingsSparsity != null)
                loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()+1);
            if(cutoffs != null)
                cutoffs.setColumnDimension(cutoffs.getColumnDimension()+1);
            if (DEBUG) {
                System.out.println("up");
            }
            increment = true;
        }
        else{
            if(loadings.getColumnDimension() == loadings.getMaxColumnDimension()) {
                from1 = -Math.log(2);
            }
            if(currentSize == 2){
                from1 = Math.log(2);
            }
            if(factors != null)
                factors.setRowDimension(factors.getRowDimension()-1);

            loadings.setColumnDimension(loadings.getColumnDimension()-1);
            if(loadingsSparsity != null)
                loadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension()-1);
            if(cutoffs != null)
                cutoffs.setColumnDimension(cutoffs.getColumnDimension()-1);
            if (DEBUG) {
                System.out.println("down");
            }
            increment = false;
        }

        //hack to let me store model state later in the code
        lfm.acceptModelState();
        if(sparsityPrior != null)
            sparsityPrior.acceptModelState();

        iterate();

        outpu = "";
//        for (int i = 0; i <storedFactors.getDimension() ; i++) {
//            outpu += storedFactors.getParameterValue(i);
//            outpu += " ";
//        }
//        System.out.println(outpu);


//        System.out.println("After");
//        System.out.println(lfm.getLogLikelihood() * (1 - sizeParam));
//        System.out.println(rowPrior.getSizeLogLikelihood());


        double finalLikelihood = lfm.getLogLikelihood() * (1 - sizeParam) + rowPrior.getSizeLogLikelihood();
//        try {
//            iterate();
//        } catch (OperatorFailedException e) {
//            e.printStackTrace();
//        }

        if (DEBUG) {
            System.out.print("storedFactors: ");
            for (int i = 0; i <storedFactors.getDimension(); i++) {
                System.out.print(storedFactors.getParameterValue(i) + " ");
            }
            System.out.println();
        }


        random = MathUtils.nextDouble();
        double test = from1 + finalLikelihood - initialLikelihood;
        test = Math.min(Math.exp(test), 1);

        boolean allRowZero = false;
//        if(loadingsSparsity != null) {
//            for (int j = 0; j < loadingsSparsity.getColumnDimension() ; j++) {
//                boolean rowZero = true;
//                for (int i = 0; i < loadingsSparsity.getRowDimension(); i++) {
//
//                    if (loadingsSparsity.getParameterValue(i, j) != 0) {
//                        rowZero = false;
//                        break;
//                    }
//
//                }
//                if(rowZero){
//                    allRowZero = true;
//                    break;
//                }
//            }
//        }


        if(random < test && (!allRowZero || !increment)){
            if (DEBUG) {
                System.out.println("accepted!\n" + test);
                System.out.println(random);
                System.out.println(test);
            }
            lfm.acceptModelState();
            lfm.makeDirty();
            if(sparsityPrior != null){
                sparsityPrior.acceptModelState();
                sparsityPrior.makeDirty();
            }
            if(loadingsPrior instanceof AbstractModelLikelihood){
                ((AbstractModelLikelihood) loadingsPrior).acceptModelState();
                loadingsPrior.makeDirty();
            }
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
            if(sparsityPrior != null)
                sparsityPrior.makeDirty();
            lfm.makeDirty();
            if(loadingsPrior instanceof AbstractModelLikelihood){
                loadingsPrior.makeDirty();
            }
            if (DEBUG) {
                System.out.println("sparsity prior lnL = " + sparsityPrior.getLogLikelihood());
                System.out.println("latentFactorModel lnL = " + lfm.getLogLikelihood());
                System.out.print("loadingsSparsity: ");
                for (int i = 0; i <loadingsSparsity.getDimension(); i++) {
                    System.out.print(loadingsSparsity.getParameterValue(i) + " ");
                }
                System.out.println();
            }
            if(factors != null)
                factors.storeParameterValues();
            loadings.storeParameterValues();
            if(loadingsSparsity != null)
                loadingsSparsity.storeParameterValues();
            if(cutoffs != null)
                cutoffs.storeParameterValues();
            lfm.acceptModelState();
            lfm.storeModelState();
            if (DEBUG) {
                System.out.println("rejected!\n" + test);
            }

        }
    }

    private void iterate() {
        if(factorOperator != null)
            factorOperator.setPathParameter(sizeParam);
        if(loadingsOperator instanceof GibbsOperator)
            loadingsOperator.setPathParameter(sizeParam);
        if(precisionGibbsOperator != null)
            precisionGibbsOperator.setPathParameter(sizeParam);
//        if(separator == null){
            separator = new double[4];
            double foWeight = 0;
            if(factorOperator != null)
                foWeight = factors.getColumnDimension() * chainLength;



            double loWeight = 0;
            if(loadingsOperator instanceof GibbsOperator){
                loWeight = loadings.getColumnDimension() * chainLength;
            }
            else{
                loWeight = loadings.getColumnDimension() * loadings.getRowDimension() * chainLength;
            }



        double sparoWeight = 0;
        if(sparsityOperator !=null)
            sparoWeight = (loadings.getRowDimension() * loadings.getColumnDimension()
//                    + factors.getColumnDimension() * factors.getRowDimension()
            ) * chainLength;
        double negWeight = 0;
        if(NOp !=null)
            negWeight = (loadings.getRowDimension() * loadings.getColumnDimension()) * chainLength;
        double precWeight = 0;
        if(precisionGibbsOperator != null){
            precWeight = chainLength;
        }

            double total = foWeight + loWeight + sparoWeight + negWeight + precWeight;
            separator[0] = foWeight / total;
            separator[1] = (foWeight + loWeight) / total;
            separator[2] = (foWeight + loWeight + sparoWeight) / total;
            separator[3] = (foWeight + loWeight + sparoWeight + negWeight)/total;
//        }
        for (int i = 0; i < total; i++) {
            double rand = MathUtils.nextDouble();
            if(rand < separator[0]){
                factorOperator.doOperation();
            }
            else if (rand < separator[1]){
                if(loadingsOperator instanceof GibbsOperator)
                    loadingsOperator.doOperation();
                else{
                    lfm.storeModelState();
                    if(loadingsPrior instanceof AbstractModelLikelihood)
                        ((AbstractModelLikelihood) loadingsPrior).storeModelState();
                    double mhRatio = - lfm.getLogLikelihood() * sizeParam - loadingsPrior.getLogLikelihood();
                    mhRatio += loadingsOperator.doOperation();
                    mhRatio += lfm.getLogLikelihood() * sizeParam + loadingsPrior.getLogLikelihood();
                    mhRatio = Math.min(1, Math.exp(mhRatio));
                    if(MathUtils.nextDouble() > mhRatio || Double.isNaN(loadingsPrior.getLogLikelihood())){
                        lfm.restoreModelState();
                        if(loadingsPrior instanceof AbstractModelLikelihood)
                            ((AbstractModelLikelihood) loadingsPrior).restoreModelState();
//                    lfm.makeDirty();
//                    sparsityPrior.makeDirty();
                    }
                    else{
                        lfm.acceptModelState();
                        if(loadingsPrior instanceof AbstractModelLikelihood)
                            ((AbstractModelLikelihood) loadingsPrior).acceptModelState();
//                    lfm.makeDirty();
//                    sparsityPrior.makeDirty();
                    }
                }
            }
            else if (rand < separator[2]){
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
            else if (rand < separator[3]){
                lfm.storeModelState();
                if (loadingsPrior instanceof AbstractModelLikelihood)
                    ((AbstractModelLikelihood) loadingsPrior).storeModelState();
                double mhRatio = -lfm.getLogLikelihood() * sizeParam - loadingsPrior.getLogLikelihood();
                mhRatio += NOp.doOperation();
                mhRatio += lfm.getLogLikelihood() * sizeParam + loadingsPrior.getLogLikelihood();
                mhRatio = Math.min(1, Math.exp(mhRatio));
                if (MathUtils.nextDouble() > mhRatio || Double.isNaN(loadingsPrior.getLogLikelihood())) {
                    lfm.restoreModelState();
                    if (loadingsPrior instanceof AbstractModelLikelihood)
                        ((AbstractModelLikelihood) loadingsPrior).restoreModelState();
//                    lfm.makeDirty();
//                    sparsityPrior.makeDirty();
                } else {
                    lfm.acceptModelState();
                    if (loadingsPrior instanceof AbstractModelLikelihood)
                        ((AbstractModelLikelihood) loadingsPrior).acceptModelState();
                }
            }
            else{
                if(precisionGibbsOperator != null)
                    precisionGibbsOperator.doOperation();
            }

        }
        if(factorOperator != null)
            factorOperator.setPathParameter(1);
        if(loadingsOperator instanceof GibbsOperator)
         loadingsOperator.setPathParameter(1);
        if(precisionGibbsOperator != null)
            precisionGibbsOperator.setPathParameter(1);
    }

    private void storeDimensions(){
        if(factors != null){
            storedFactors.setRowDimension(factors.getRowDimension());
            storedFactors.setColumnDimension(factors.getColumnDimension());
        }
        storedLoadings.setRowDimension(loadings.getRowDimension());
        storedLoadings.setColumnDimension(loadings.getColumnDimension());
        if(loadingsSparsity != null){
            storedLoadingsSparsity.setRowDimension(loadingsSparsity.getRowDimension());
            storedLoadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension());
        }
        if(cutoffs != null){
            storedCutoffs.setRowDimension(cutoffs.getRowDimension());
            storedCutoffs.setColumnDimension(cutoffs.getColumnDimension());
        }
    }

    private void restoreDimensions(){
        if(factors != null){
            factors.setRowDimension(storedFactors.getRowDimension());
            factors.setColumnDimension(storedFactors.getColumnDimension());
        }
        loadings.setRowDimension(storedLoadings.getRowDimension());
        loadings.setColumnDimension(storedLoadings.getColumnDimension());
        if(loadingsSparsity != null){
            loadingsSparsity.setRowDimension(storedLoadingsSparsity.getRowDimension());
            loadingsSparsity.setColumnDimension(storedLoadingsSparsity.getColumnDimension());
        }
        if(cutoffs != null){
            cutoffs.setRowDimension(storedCutoffs.getRowDimension());
            cutoffs.setColumnDimension(storedCutoffs.getColumnDimension());
        }
    }


    private void storeValues(){
        if(factors != null){
        for (int i = 0; i < factors.getDimension(); i++) {
            storedFactors.setParameterValue(i, factors.getParameterValue(i));}
        }

        for (int i = 0; i < loadings.getDimension(); i++) {
            storedLoadings.setParameterValue(i, loadings.getParameterValue(i));
            if(loadingsSparsity != null)
                storedLoadingsSparsity.setParameterValue(i, loadingsSparsity.getParameterValue(i));
            if(storedCutoffs != null)
               storedCutoffs.setParameterValue(i, cutoffs.getParameterValue(i));
        }
    }

    private void restoreValues(){
        if(factors != null){
        for (int i = 0; i < factors.getDimension(); i++) {
            factors.setParameterValue(i, storedFactors.getParameterValue(i));}
        }
        for (int i = 0; i < loadings.getDimension(); i++) {
            loadings.setParameterValue(i, storedLoadings.getParameterValue(i));
            if(loadingsSparsity != null)
                loadingsSparsity.setParameterValue(i, storedLoadingsSparsity.getParameterValue(i));
            if(cutoffs != null)
                cutoffs.setParameterValue(i, storedCutoffs.getParameterValue(i));
        }
    }
}
