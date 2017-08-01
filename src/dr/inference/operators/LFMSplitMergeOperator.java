package dr.inference.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NormalDistribution;

/**
 * @author Max Tolkoff
 */
public class LFMSplitMergeOperator extends AbstractCoercableOperator {

    NormalDistribution drawDistribution;
    AdaptableSizeFastMatrixParameter factors;
    AdaptableSizeFastMatrixParameter sparseLoadings;
    AdaptableSizeFastMatrixParameter denseLoadings;
    AdaptableSizeFastMatrixParameter cutoffs;
    FullyConjugateMultivariateTraitLikelihood tree;
    GaussianProcessFromTree treeDraw;
    NormalDistribution standardNormal = new NormalDistribution(0, 1);
    GammaDistribution gamma = new GammaDistribution(1, 1);

    private String lastCall;

    //TODO implment flag for not this case
    final boolean upperTriangular = true;

    public LFMSplitMergeOperator(double weight, double initialVariance, AdaptableSizeFastMatrixParameter factors,
                                 AdaptableSizeFastMatrixParameter denseLoadings, AdaptableSizeFastMatrixParameter sparseLoadings,
                                 AdaptableSizeFastMatrixParameter cutoffs, FullyConjugateMultivariateTraitLikelihood tree){

        super(CoercionMode.DEFAULT);
        setWeight(weight);

        drawDistribution = new NormalDistribution(0, Math.sqrt(initialVariance));
        this.factors = factors;
        this.sparseLoadings = sparseLoadings;
        this.denseLoadings = denseLoadings;
        this.cutoffs = cutoffs;
        this.tree = tree;
        GaussianProcessFromTree treeDraw = new GaussianProcessFromTree(tree);
        this.treeDraw = treeDraw;
    }


    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Model Split-Merge Operator";
    }

    @Override
    public double doOperation() {
        double randNum = MathUtils.nextDouble();
        double from1 = 0;
        double hastingsRatio;
        if((randNum < .5 ||
                denseLoadings.getColumnDimension() == denseLoadings.getMaxColumnDimension()) && denseLoadings.getColumnDimension() != 1){
            lastCall = "decrement";
            if(
                    denseLoadings.getColumnDimension() == 2
                    )
                from1 += Math.log(2);
            if(
                    denseLoadings.getColumnDimension() == denseLoadings.getMaxColumnDimension()
                    )
                from1 += -Math.log(2);
            hastingsRatio = decrement();
//            hastingsRatio = Double.POSITIVE_INFINITY;
        }
        else{
            lastCall = "increment";
            if(
                    denseLoadings.getColumnDimension() == 1
                    )
                from1 += -Math.log(2);
            if(
                    denseLoadings.getColumnDimension() == denseLoadings.getMaxColumnDimension() - 1
                    )
                from1 += Math.log(2);
            hastingsRatio = increment();
//            hastingsRatio = Double.POSITIVE_INFINITY;
        }

        denseLoadings.fireParameterChangedEvent();
        sparseLoadings.fireParameterChangedEvent();
        cutoffs.fireParameterChangedEvent();
        factors.fireParameterChangedEvent();
//        factors.fireParameterChangedEvent(-1, Variable.ChangeType.ADDED);

//        System.out.println(drawDistribution.getSD());
//        System.out.println(hastingsRatio);
        return hastingsRatio + from1;
    }

    @Override
    public double getCoercableParameter() {
//        return Math.log(drawDistribution.getSD());
        return 0;
    }

    @Override
    public void setCoercableParameter(double value) {
//        drawDistribution.setSD(Math.exp(value));

    }

    @Override
    public double getRawParameter() {
//        return drawDistribution.getSD();
        return 0;
    }

    private double increment(){
        double hastings = 0;
        int splitFac = MathUtils.nextInt(denseLoadings.getColumnDimension());
        int lastColumn = denseLoadings.getColumnDimension();
        denseLoadings.setColumnDimension(denseLoadings.getColumnDimension() + 1);
        sparseLoadings.setColumnDimension(sparseLoadings.getColumnDimension() + 1);
        cutoffs.setColumnDimension(cutoffs.getColumnDimension() + 1);
//        factors.setRowDimension(factors.getRowDimension() + 1);
//        System.out.println(factors.getColumnDimension());
//        NormalDistribution standardNormal = new NormalDistribution(0, 1);
//        double draw = 0;
//        double cutoffDraw = 0;

        for (int i = 0; i < denseLoadings.getRowDimension(); i++) {
            if(i >= lastColumn || !upperTriangular) {
//                hastings += Math.log(3);

                if (sparseLoadings.getParameterValue(i, splitFac) == 0) {
                    hastings += sparseIncrement(i, lastColumn);
                } else {
                    //discrete split
                    int intDraw = MathUtils.nextInt(3);
                    intDraw = 2;
                    if (intDraw == 0) {
                        hastings += sparseIncrement(i, lastColumn);
                    } else if (intDraw == 1) {
                        //move to last column
                        sparseLoadings.setParameterValueQuietly(i, lastColumn, sparseLoadings.getParameterValue(i, splitFac));
                        cutoffs.setParameterValueQuietly(i, lastColumn, cutoffs.getParameterValue(i, splitFac));
                        denseLoadings.setParameterValueQuietly(i, lastColumn, denseLoadings.getParameterValue(i, splitFac));

                        //split
                        hastings += sparseIncrement(i, splitFac);
                    } else {

                        //backwards cutoff draw hastings ratio
                        hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, splitFac), 2));
                        //both are 1
                        sparseLoadings.setParameterValueQuietly(i, lastColumn, 1);

                        //dense split
                        double draw = (Double) drawDistribution.nextRandom();
                        double reset1 = denseLoadings.getParameterValue(i, splitFac) + draw;
                        double reset2 = denseLoadings.getParameterValue(i, splitFac) - draw;
                        denseLoadings.setParameterValueQuietly(i, lastColumn, reset1);
                        denseLoadings.setParameterValueQuietly(i, splitFac, reset2);
//                        System.out.println(draw);
                        hastings += -drawDistribution.logPdf(draw);
                        //jacobian
                        hastings += Math.log(2);

                        //new cutoffs
                        double cutoffDraw1 = MathUtils.nextDouble() * Math.pow(reset1, 2);
                        double cutoffDraw2 = MathUtils.nextDouble() * Math.pow(reset2, 2);
                        hastings += Math.log(Math.pow(reset1, 2));
                        hastings += Math.log(Math.pow(reset2, 2));
                        cutoffs.setParameterValueQuietly(i, lastColumn, cutoffDraw1);
                        cutoffs.setParameterValueQuietly(i, splitFac, cutoffDraw2);
                    }
                }
            }
            else{
//                hastings += -standardNormal.logPdf(0);
                double cutoffDraw = gamma.nextGamma();
                hastings += -gamma.logPdf(cutoffDraw);
                cutoffs.setParameterValueQuietly(i, lastColumn, cutoffDraw);
//                sparseLoadings.setParameterValueQuietly(i, lastColumn, MathUtils.nextInt(2));
//                hastings += Math.log(2);
            }

        }
        //setup tree draw
//        double[] drawnFactors = treeDraw.nextRandomFast();
////        System.out.println(factors.getColumnDimension());
////        System.out.println(drawnFactors.length);
//        for (int i = 0; i < factors.getColumnDimension() ; i++) {
////            System.out.println(i);
//            factors.setParameterValueQuietly(lastColumn, i, drawnFactors[i]);
//        }
//        factors.fireParameterChangedEvent();
//        double logFull = tree.getLogLikelihood();
//        factors.setRowDimension(factors.getRowDimension() - 1);
//        double logEmpty = tree.getLogLikelihood();
//        hastings += logEmpty - logFull;
////        System.out.println(logEmpty - logFull);
//        factors.setRowDimension(factors.getRowDimension() + 1);
//
//        for (int i = 0; i < factors.getColumnDimension(); i++) {
////            double draw = (Double) drawDistribution.nextRandom();
//            factors.setParameterValueQuietly(lastColumn, i, factors.getParameterValue(splitFac, i)  + drawnFactors[i] * drawDistribution.getSD());
//            factors.setParameterValueQuietly(splitFac, i, factors.getParameterValue(splitFac, i)   - drawnFactors[i] * drawDistribution.getSD());
////            hastings += -drawDistribution.logPdf(draw);
//            //Jacobian
//            hastings += Math.log(2);
//        }

//        System.out.println(factors.getRowDimension());
        return hastings;
    }

    double sparseIncrement(int i, int column){

        double draw = 0;
        double cutoffDraw = 0;
        double hastings = 0;

        //discrete split
        sparseLoadings.setParameterValueQuietly(i, column, 0);
//        hastings += Math.log(3);

        //dense denovo draw
        draw = (Double) standardNormal.nextRandom();
        denseLoadings.setParameterValueQuietly(i, column, draw);
        hastings += -standardNormal.logPdf(draw);

        //cutoff draw
        cutoffDraw = MathUtils.nextDouble() * Math.pow(draw, 2);
        hastings += Math.log(Math.pow(draw, 2));
        cutoffs.setParameterValueQuietly(i, column, cutoffDraw);
//        System.out.println(draw);
        return hastings;
    }


    private double decrement(){
        double hastings = 0;
        int combineFac = MathUtils.nextInt(denseLoadings.getColumnDimension() - 1);
        int lastColumn = denseLoadings.getColumnDimension() - 1;
        for (int i = 0; i < denseLoadings.getRowDimension(); i++) {
            if(i >= lastColumn || !upperTriangular){
//                System.out.println(i);
                //Backwards probabilities
//                hastings += -Math.log(3);
                if(sparseLoadings.getParameterValue(i, combineFac) == 1 && sparseLoadings.getParameterValue(i, lastColumn) == 1){
                    hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, combineFac), 2));
                    hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, lastColumn), 2));
                    hastings += drawDistribution.logPdf(( denseLoadings.getParameterValue(i, lastColumn) - denseLoadings.getParameterValue(i, combineFac)) / 2);
                    //jacobian
                    hastings += -Math.log(2);
                    double sum = (denseLoadings.getParameterValue(i, combineFac) + denseLoadings.getParameterValue(i, lastColumn)) / 2;
                    denseLoadings.setParameterValueQuietly(i, combineFac, sum);
                    double draw = MathUtils.nextDouble() * Math.pow(sum, 2);
                    cutoffs.setParameterValueQuietly(i, combineFac, draw);
                }
                else if(sparseLoadings.getParameterValue(i, combineFac) == 1 && sparseLoadings.getParameterValue(i, lastColumn) == 0){
                    hastings += standardNormal.logPdf(denseLoadings.getParameterValue(i, lastColumn));
                    hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, lastColumn), 2));
                }
                else if(sparseLoadings.getParameterValue(i, combineFac) == 0 && sparseLoadings.getParameterValue(i, lastColumn) == 1){
                    hastings += standardNormal.logPdf(denseLoadings.getParameterValue(i, combineFac));
                    hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, combineFac), 2));
                    denseLoadings.setParameterValueQuietly(i, combineFac, denseLoadings.getParameterValue(i, lastColumn));
                    cutoffs.setParameterValueQuietly(i, combineFac, cutoffs.getParameterValue(i, lastColumn));
                    sparseLoadings.setParameterValueQuietly(i, combineFac, 1);
                }
                else{
                    hastings += standardNormal.logPdf(denseLoadings.getParameterValue(i, lastColumn));
                    hastings += -Math.log(Math.pow(denseLoadings.getParameterValue(i, lastColumn), 2));
                }
            }
            else{
//                hastings += standardNormal.logPdf(0);
                hastings += gamma.logPdf(cutoffs.getParameterValue(i, lastColumn));
//                hastings += -Math.log(2);
            }
        }



//        for (int i = 0; i < factors.getColumnDimension(); i++) {
////            hastings += drawDistribution.logPdf((factors.getParameterValue(combineFac, i) - factors.getParameterValue(lastColumn, i)) / 2 );
//            hastings += -Math.log(2);
////            System.out.println(drawDistribution.logPdf((factors.getParameterValue(combineFac, i) - factors.getParameterValue(lastColumn, i)) / 2 ));
//            double temp1 = (factors.getParameterValue(combineFac, i) + factors.getParameterValue(lastColumn, i)) / 2 ;
//            double temp2 = (factors.getParameterValue(combineFac, i) - factors.getParameterValue(lastColumn, i)) / (2 * drawDistribution.getSD());
//            factors.setParameterValueQuietly(combineFac, i, temp1 );
//            factors.setParameterValueQuietly(lastColumn, i, temp2 );
//        }
//        factors.fireParameterChangedEvent();
//        double fullDraw = tree.getLogLikelihood();
//
//
        denseLoadings.setColumnDimension(denseLoadings.getColumnDimension() - 1);
        sparseLoadings.setColumnDimension(sparseLoadings.getColumnDimension() - 1);
        cutoffs.setColumnDimension(cutoffs.getColumnDimension() - 1);
////        System.out.println(factors.getRowDimension());
//        factors.setRowDimension(factors.getRowDimension() - 1);
//////
//        double truncDraw = tree.getLogLikelihood();
//        hastings += fullDraw - truncDraw;


//        System.out.println(fullDraw - truncDraw);
//        System.out.println(hastings);
        return hastings;
    }


//    @Override
//    public void accept(double deviation) {
//        super.accept(deviation);
//        System.out.println(lastCall);
//        System.out.println(deviation);
//    }

//    public void reject(){
//        super.reject();
//        System.out.println(lastCall);
//    }

    //    private void storeDimensions(){
//        if(factors != null){
//            storedFactors.setRowDimension(factors.getRowDimension());
//            storedFactors.setColumnDimension(factors.getColumnDimension());
//        }
//        storedLoadings.setRowDimension(loadings.getRowDimension());
//        storedLoadings.setColumnDimension(loadings.getColumnDimension());
//        if(loadingsSparsity != null){
//            storedLoadingsSparsity.setRowDimension(loadingsSparsity.getRowDimension());
//            storedLoadingsSparsity.setColumnDimension(loadingsSparsity.getColumnDimension());
//        }
//        if(cutoffs != null){
//            storedCutoffs.setRowDimension(cutoffs.getRowDimension());
//            storedCutoffs.setColumnDimension(cutoffs.getColumnDimension());
//        }
//    }
//
//    private void restoreDimensions(){
//        if(factors != null){
//            factors.setRowDimension(storedFactors.getRowDimension());
//            factors.setColumnDimension(storedFactors.getColumnDimension());
//        }
//        loadings.setRowDimension(storedLoadings.getRowDimension());
//        loadings.setColumnDimension(storedLoadings.getColumnDimension());
//        if(loadingsSparsity != null){
//            loadingsSparsity.setRowDimension(storedLoadingsSparsity.getRowDimension());
//            loadingsSparsity.setColumnDimension(storedLoadingsSparsity.getColumnDimension());
//        }
//        if(cutoffs != null){
//            cutoffs.setRowDimension(storedCutoffs.getRowDimension());
//            cutoffs.setColumnDimension(storedCutoffs.getColumnDimension());
//        }
//    }
//
//
//    private void storeValues(){
//        if(factors != null){
//            for (int i = 0; i < factors.getDimension(); i++) {
//                storedFactors.setParameterValue(i, factors.getParameterValue(i));}
//        }
//
//        for (int i = 0; i < loadings.getDimension(); i++) {
//            storedLoadings.setParameterValue(i, loadings.getParameterValue(i));
//            if(loadingsSparsity != null)
//                storedLoadingsSparsity.setParameterValue(i, loadingsSparsity.getParameterValue(i));
//            if(storedCutoffs != null)
//                storedCutoffs.setParameterValue(i, cutoffs.getParameterValue(i));
//        }
//    }
//
//    private void restoreValues(){
//        if(factors != null){
//            for (int i = 0; i < factors.getDimension(); i++) {
//                factors.setParameterValue(i, storedFactors.getParameterValue(i));}
//        }
//        for (int i = 0; i < loadings.getDimension(); i++) {
//            loadings.setParameterValue(i, storedLoadings.getParameterValue(i));
//            if(loadingsSparsity != null)
//                loadingsSparsity.setParameterValue(i, storedLoadingsSparsity.getParameterValue(i));
//            if(cutoffs != null)
//                cutoffs.setParameterValue(i, storedCutoffs.getParameterValue(i));
//        }
//    }


}
