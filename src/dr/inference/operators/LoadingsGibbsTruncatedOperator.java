package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import jebl.math.Random;

import java.util.ArrayList;
import java.util.ListIterator;

/**
 * Created by max on 2/4/16.
 */
public class LoadingsGibbsTruncatedOperator extends SimpleMCMCOperator implements GibbsOperator{
    MomentDistributionModel prior;
    LatentFactorModel LFM;
    double[][] precisionArray;
    double[] meanMidArray;
    double[] meanArray;
    boolean randomScan;
    double pathParameter=1.0;


    double priorPrecision;
    double priorMeanPrecision;
    MatrixParameterInterface loadings;

    DistributionLikelihood cutoffPrior;

    public LoadingsGibbsTruncatedOperator(LatentFactorModel LFM, MomentDistributionModel prior, double weight, boolean randomScan, MatrixParameterInterface loadings, DistributionLikelihood cutoffPrior) {
        setWeight(weight);

        this.loadings=loadings;
        this.prior = prior;
        this.LFM = LFM;
        priorPrecision = (this.prior.getScaleMatrix()[0][0]);
        priorMeanPrecision = this.prior.getMean()[0] * priorPrecision;
        this.cutoffPrior = cutoffPrior;
    }

    private void getPrecisionOfTruncated(MatrixParameterInterface full, int newRowDimension, int row, double[][] answer) {
        int p = full.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += full.getParameterValue(i, k) * full.getParameterValue(j, k);
                answer[i][j] = sum * LFM.getColumnPrecision().getParameterValue(row, row);
                if (i == j) {
                    answer[i][j] =answer[i][j] * pathParameter + priorPrecision;
                } else {
                    answer[i][j] *= pathParameter;
                    answer[j][i] = answer[i][j];
                }
            }
        }
    }


    private void getTruncatedMean(int newRowDimension, int dataColumn, double[][] variance, double[] midMean, double[] mean) {
        MatrixParameterInterface data = LFM.getScaledData();
        MatrixParameterInterface Left = LFM.getFactors();
        int p = data.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
            for (int k = 0; k < p; k++)
                sum += Left.getParameterValue(i, k) * data.getParameterValue(dataColumn, k);
            sum = sum * LFM.getColumnPrecision().getParameterValue(dataColumn, dataColumn);
            sum += priorMeanPrecision;
            midMean[i] = sum;
        }
        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
            for (int k = 0; k < newRowDimension; k++)
                sum += variance[i][k] * midMean[k];
            mean[i] = sum;
        }

    }

    private void getPrecision(int i, double[][] answer) {
        int size = LFM.getFactorDimension();
        getPrecisionOfTruncated(LFM.getFactors(), size, i, answer);
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean) {
        int size = LFM.getFactorDimension();
        getTruncatedMean(size, i, variance, midMean, mean);
        for (int j = 0; j <mean.length ; j++) {//TODO implement for generic prior
            mean[j]*=pathParameter;
        }

    }

    private void copy(int i, double[] random) {
        MatrixParameterInterface changing = loadings;
        for (int j = 0; j < random.length; j++) {
            changing.setParameterValueQuietly(i, j, random[j]);
        }
    }

    private void getTruncatedDraw(int row, int column, NormalDistribution truncated){
        double lowCutoff;
        double highCutoff;
        MatrixParameterInterface cutoff = (MatrixParameterInterface) prior.getCutoff();
        lowCutoff = - Math.sqrt(cutoff.getParameterValue(row, column));
        highCutoff = - lowCutoff;
        double low = truncated.cdf(lowCutoff);
        double high = truncated.cdf(highCutoff);
        double split = low / (low + (1-high));
        double rand = MathUtils.nextDouble();
        double draw;
        if(rand < split){
            draw = MathUtils.nextDouble() * low;
            draw = truncated.quantile(draw);
        }
        else{
            draw = MathUtils.nextDouble() * (1- high) + high;
            draw = truncated.quantile(draw);
        }


        if(!Double.isNaN(draw)){
            loadings.setParameterValue(row, column, draw);
        }

    }

    private void drawI(int i, int column) {
        double[] draws = null;
        precisionArray = new double[loadings.getColumnDimension()][loadings.getColumnDimension()];
        double[][] variance;
        meanMidArray = new double[loadings.getColumnDimension()];
        meanArray = new double[loadings.getColumnDimension()];
        double[][] cholesky = null;
        NormalDistribution conditioned;
        getPrecision(i, precisionArray);
        if(LFM.getLoadings().getParameterValue(i, column) != 0){
        variance = (new SymmetricMatrix(precisionArray)).inverse().toComponents();


//        try {
//            cholesky = new CholeskyDecomposition(variance).getL();
//        } catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();
//        }

        getMean(i, variance, meanMidArray, meanArray);


        conditioned = getConditionalDistribution(meanArray, variance, column, i);
        }
        else
            conditioned = new NormalDistribution(0, Math.sqrt(1 / prior.getScaleMatrix()[0][0]));


        if(MathUtils.nextDouble() < .5) {
            getTruncatedDraw(i, column, conditioned);
            getCutoffDraw(i, column);
        }
        else{
            getCutoffDraw(i, column);
            getTruncatedDraw(i, column, conditioned);
        }
    }

    private NormalDistribution getConditionalDistribution(double[] meanArray, double[][] variance, int column, int row) {
        double[][] newVariance = new double[meanArray.length - 1][meanArray.length - 1];
        for (int i = 0; i < meanArray.length; i++) {
            for (int j = 0; j < meanArray.length; j++) {
                if(i < column && j < column){
                    newVariance[i][j] = variance[i][j];
                }
                else if(i < column && j > column){
                    newVariance[i][j - 1] = variance[i][j];
                }
                else if(i > column && j < column){
                    newVariance[i - 1][j] = variance[i][j];
                }
                else if(i > column && j > column){
                    newVariance[i - 1][j - 1] = variance[i][j];
                }
                else{}
            }
        }
        double[][] smallPrecision = (new SymmetricMatrix(newVariance)).inverse().toComponents();
        double[] meanStore1 = new double[meanArray.length - 1];
        double[] meanStore2 = new double[meanArray.length - 1];
        double[] precStore = new double[meanArray.length - 1];
        for (int i = 0; i < meanArray.length; i++) {
            if(i < column){
                meanStore1[i] = LFM.getLoadings().getParameterValue(row, i) - meanArray[i];
            }
            else if (i > column){
                meanStore1[i - 1] = LFM.getLoadings().getParameterValue(row, i) - meanArray[i];
            }
            else{}
        }
        for (int i = 0; i < meanArray.length - 1; i++) {
            for (int j = 0; j < meanArray.length - 1; j++) {
                meanStore2[i] += smallPrecision[i][j] * meanStore1[j];
            }
        }
        double mean = meanArray[column];
        for (int i = 0; i < meanArray.length - 1; i++) {
            if(i < column){
                mean += meanStore2[i] * variance[i][column];
            }
            else{
                mean += meanStore2[i] * variance[i + 1][column];
            }
        }
        for (int i = 0; i < meanArray.length - 1; i++) {
            for (int j = 0; j < meanArray.length - 1; j++) {
                if(i < column)
                    precStore[i] += smallPrecision[i][j] * variance[j][column];
                else
                    precStore[i] += smallPrecision[i][j] * variance[j+1][column];
            }

        }
        double varianceElement = variance[column][column];
        for (int i = 0; i < meanArray.length - 1; i++) {
            if(i < column)
                varianceElement -= precStore[i] * variance[i][column];
            else
                varianceElement -= precStore[i] * variance[i+1][column];
        }

        return new NormalDistribution(mean, Math.sqrt(varianceElement));



    }

    void getCutoffDraw(int row, int column){
        double loadingsCutoff = Math.abs(loadings.getParameterValue(row, column));
        double stopperCDF = Math.pow(cutoffPrior.getDistribution().cdf(loadingsCutoff), 2);
        double randQuant = MathUtils.nextDouble() * stopperCDF;
        ((MatrixParameterInterface) prior.getCutoff()).setParameterValue(row, column, cutoffPrior.getDistribution().quantile(randQuant));
    }

    @Override
    public int getStepCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return "loadingsGibbsTruncatedOperator";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {

        int size = LFM.getLoadings().getRowDimension();
        int column = MathUtils.nextInt(LFM.getLoadings().getColumnDimension());
            for (int i = 0; i < size; i++) {
                drawI(i, column);
            }
            ((Parameter) loadings).fireParameterChangedEvent();
        return 0;
    }

    public void setPathParameter(double beta){
        pathParameter=beta;
    }
}

