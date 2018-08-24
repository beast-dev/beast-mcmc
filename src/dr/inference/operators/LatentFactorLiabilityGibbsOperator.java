package dr.inference.operators;

import dr.evomodel.continuous.OrderedLatentLiabilityLikelihood;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

/**
 * Created by Max on 9/1/16.
 */
public class LatentFactorLiabilityGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{

    LatentFactorModel lfm;
    OrderedLatentLiabilityLikelihood liabilityLikelihood;

    public LatentFactorLiabilityGibbsOperator(double weight, LatentFactorModel lfm, OrderedLatentLiabilityLikelihood liabilityLikelihood) {
        setWeight(weight);


        this.lfm = lfm;
        this.liabilityLikelihood = liabilityLikelihood;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "LatentFactorLiabilityGibbsOperator";
    }

    @Override
    public double doOperation() {
        if(liabilityLikelihood.getOrdering())
            doUnorderedOperation();
        else
            doOrderedOperation();

        return 0;
    }

    void doUnorderedOperation() {
        double[] LxF = lfm.getLxF();
        DiagonalMatrix colPrec = (DiagonalMatrix) lfm.getColumnPrecision();
        Parameter continuous = lfm.getContinuous();
        MatrixParameterInterface lfmData = lfm.getScaledData();

        for (int i = 0; i < lfmData.getColumnDimension(); i++) {
            int LLpointer = 0;
            int[] data = liabilityLikelihood.getData(i);
            for (int index = 0; index < data.length; ++index) {

                int datum = data[index];
                Parameter numClasses = liabilityLikelihood.numClasses;
                int dim = (int) numClasses.getParameterValue(index);
                if(datum >= dim  && continuous.getParameterValue(LLpointer) == 0){
                    double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + LLpointer], colPrec.getParameterValue(LLpointer, LLpointer), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                    lfmData.setParameterValue(LLpointer, i, draw);
                }

                if (dim == 1.0) {
                    if (continuous.getParameterValue(LLpointer) == 0) {
                        double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + LLpointer], colPrec.getParameterValue(LLpointer, LLpointer), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                        lfmData.setParameterValue(LLpointer, i, draw);
                    }

                    LLpointer++;
                } else if (dim == 2.0) {


                    if (datum == 0) {
                        double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + LLpointer], colPrec.getParameterValue(LLpointer, LLpointer), Double.NEGATIVE_INFINITY, 0);
                        lfmData.setParameterValue(LLpointer, i, draw);

                    } else {
                        double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + LLpointer], colPrec.getParameterValue(LLpointer, LLpointer), 0, Double.POSITIVE_INFINITY);
                        lfmData.setParameterValue(LLpointer, i, draw);

                    }
                    LLpointer++;
                } else {
                    double[] trait = new double[dim];
                    trait[0] = 0.0;
                    if(datum == 0){
                        for (int l = 0; l < dim - 1; l++) {
                             lfmData.setParameterValue(LLpointer + l, i, drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + (LLpointer + l)], colPrec.getParameterValue((LLpointer + l), (LLpointer + l)), Double.NEGATIVE_INFINITY, 0));
                        }
                    }
                    else {
                        trait[datum] = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + (LLpointer + datum - 1)], colPrec.getParameterValue((LLpointer + datum - 1), (LLpointer + datum - 1)), 0, Double.POSITIVE_INFINITY);
                        lfmData.setParameterValue(LLpointer + datum - 1, i, trait[datum]);
                        for (int l = 1; l < dim; l++) {
                            if(l != datum){
//                                System.out.println("Free Rolls");
//                                System.out.println(LxF[i * lfmData.getRowDimension() + (LLpointer + l - 1)]);
                                trait[l] = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + (LLpointer + l - 1)], colPrec.getParameterValue((LLpointer + l - 1), (LLpointer + l - 1)), Double.NEGATIVE_INFINITY, trait[datum]);
                                lfmData.setParameterValue(LLpointer + l - 1, i, trait[l]);
                            }
                        }
//                        double max = Double.NEGATIVE_INFINITY;
//                        for (int j = 0; j < trait.length; j++) {
//                            if(max < trait[j]){
//                                max = trait[j];
//                            }
//                        }
//                        System.out.println("Constrained");
//                        System.out.println(LxF[i * lfmData.getRowDimension() + (LLpointer + datum - 1)]);


                    }
//                    valid = isMax(trait, datum);


                    LLpointer += dim - 1;

                }
            }
        }
    }



    void doOrderedOperation(){
        double[] LxF = lfm.getLxF();
        DiagonalMatrix colPrec = (DiagonalMatrix) lfm.getColumnPrecision();
        Parameter continuous = lfm.getContinuous();
        MatrixParameterInterface lfmData = lfm.getScaledData();
        Parameter threshold = liabilityLikelihood.getThreshold();

        for (int i = 0; i < lfmData.getColumnDimension(); i++) {
            int Thresholdpointer = 0;
            int[] data = liabilityLikelihood.getData(i);
            for (int index = 0; index < data.length; ++index) {

                int datum = data[index];
                Parameter numClasses = liabilityLikelihood.numClasses;
                int dim = (int) numClasses.getParameterValue(index);
                if(datum >= dim  && continuous.getParameterValue(index) == 0){
                    double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + index], colPrec.getParameterValue(index, index), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                    lfmData.setParameterValue(index, i, draw);
                }
                else {
                    if (dim == 1.0) {
                        if (continuous.getParameterValue(index) == 0) {
                            double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + index], colPrec.getParameterValue(index, index), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
                            lfmData.setParameterValue(index, i, draw);
                        }
                    } else if (dim == 2.0) {


                        if (datum == 0) {
                            double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + index], colPrec.getParameterValue(index, index), Double.NEGATIVE_INFINITY, 0);
                            lfmData.setParameterValue(index, i, draw);

                        } else {
                            double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + index], colPrec.getParameterValue(index, index), 0, Double.POSITIVE_INFINITY);
                            lfmData.setParameterValue(index, i, draw);

                        }
                    } else {
                        double[] thresholdList = new double[dim + 1];
                        thresholdList[0] = Double.NEGATIVE_INFINITY;
                        thresholdList[1] = 0;
                        thresholdList[dim] = Double.POSITIVE_INFINITY;

                        for (int j = 0; j < thresholdList.length - 3; j++) {
                            thresholdList[j + 2] = threshold.getParameterValue(Thresholdpointer + j);
                        }
                        Thresholdpointer += dim - 2;
                        double draw = drawTruncatedNormalDistribution(LxF[i * lfmData.getRowDimension() + index], colPrec.getParameterValue(index, index), thresholdList[datum], thresholdList[datum + 1]);
                        lfmData.setParameterValue(index, i, draw);

                    }
                }
//                    valid = isMax(trait, datum);



            }
        }
    }




    double drawTruncatedNormalDistribution(double mean, double precision, double lower, double upper){
        double sd = Math.sqrt(1 / precision);
        NormalDistribution normal = new NormalDistribution(mean, sd);
        double newLower = normal.cdf(lower);
        double newUpper = normal.cdf(upper);


        double cdfDraw = 1.0;
        int iterator = 0;
        boolean invalid = true;
        double draw = 0;
        while(iterator < 10000 && invalid){
                cdfDraw = MathUtils.nextDouble() * (newUpper - newLower) + newLower;
                draw = normal.quantile(cdfDraw);
                if(!Double.isNaN(draw) && draw > lower && draw < upper) {
                    invalid = false;
                }
                iterator++;
        }
//        if(iterator != 1){
//            System.out.println(iterator);
//            System.out.println(draw);
//            System.out.println(lower);
//            System.out.println(upper);}
        if(Double.isNaN(draw) || Double.isInfinite(draw)){
            if(Double.isInfinite(lower)){
//                System.out.println("upper");
//                System.out.println(upper);
                return upper;}
            else if(Double.isInfinite(upper)){
//                System.out.println("lower");
//                System.out.println(lower);
                return lower;}
            else
                return (lower + upper) / 2;
        }
        else
            return draw;
    }
}
