package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MomentDistributionModel;
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
    ArrayList<double[][]> precisionArray;
    ArrayList<double[]> meanMidArray;
    ArrayList<double[]> meanArray;
    boolean randomScan;
    double pathParameter=1.0;


    double priorPrecision;
    double priorMeanPrecision;
    MatrixParameterInterface loadings;

    public LoadingsGibbsTruncatedOperator(LatentFactorModel LFM, MomentDistributionModel prior, double weight, boolean randomScan, MatrixParameterInterface loadings) {
        setWeight(weight);

        this.loadings=loadings;
        this.prior = prior;
        this.LFM = LFM;
        precisionArray = new ArrayList<double[][]>();
        double[][] temp;
        this.randomScan = randomScan;


        meanArray = new ArrayList<double[]>();
        meanMidArray = new ArrayList<double[]>();
        double[] tempMean;
        if (!randomScan) {
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                temp = new double[i + 1][i + 1];
                precisionArray.add(temp);
            }
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[i + 1];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[i + 1];
                meanMidArray.add(tempMean);
            }
        } else {
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                temp = new double[LFM.getFactorDimension() - i][LFM.getFactorDimension() - i];
                precisionArray.add(temp);
            }
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[LFM.getFactorDimension() - i];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[LFM.getFactorDimension() - i];
                meanMidArray.add(tempMean);
            }
        }

//            vectorProductAnswer=new MatrixParameter[LFM.getLoadings().getRowDimension()];
//            for (int i = 0; i <vectorProductAnswer.length ; i++) {
//                vectorProductAnswer[i]=new MatrixParameter(null);
//                vectorProductAnswer[i].setDimensions(i+1, 1);
//            }

//        priorMeanVector=new MatrixParameter[LFM.getLoadings().getRowDimension()];
//            for (int i = 0; i <priorMeanVector.length ; i++) {
//                priorMeanVector[i]=new MatrixParameter(null, i+1, 1, this.prior.getMean()/(this.prior.getSD()*this.prior.getSD()));
//
//
//            }
        priorPrecision = (this.prior.getScaleMatrix()[0][0]);
        priorMeanPrecision = this.prior.getMean()[0] * priorPrecision;
    }

    private void getPrecisionOfTruncated(MatrixParameterInterface full, int newRowDimension, int row, double[][] answer) {

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        int p = full.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += full.getParameterValue(i, k) * full.getParameterValue(j, k);
                answer[i][j] = sum * LFM.getColumnPrecision().getParameterValue(row, row);
                if (i == j) {
                    answer[i][j] =answer[i][j]*pathParameter+ priorPrecision;
                } else {
                    answer[i][j]*=pathParameter;
                    answer[j][i] = answer[i][j];
                }
            }
        }
    }


    private void getTruncatedMean(int newRowDimension, int dataColumn, double[][] variance, double[] midMean, double[] mean) {

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());
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
        if (i < size) {
            getPrecisionOfTruncated(LFM.getFactors(), i + 1, i, answer);
        } else {
            getPrecisionOfTruncated(LFM.getFactors(), size, i, answer);
        }
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean) {
//        Matrix factors=null;
        int size = LFM.getFactorDimension();
//        double[] scaledDataColumn=LFM.getScaledData().getRowValues(i);
//        Vector dataColumn=null;
//        Vector priorVector=null;
//        Vector temp=null;
//        Matrix data=new Matrix(LFM.getScaledData().getParameterAsMatrix());
        if (i < size) {
            getTruncatedMean(i + 1, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[i].add(vectorProductAnswer[i]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        } else {
            getTruncatedMean(size, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[size-1].add(vectorProductAnswer[size-1]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        }
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

    private double[] getDraws(int row, double[] mean, double[][] Cholesky){
        double[] temp = new double[mean.length];
        double[] draws = new double[mean.length];
        double lowCutoff;
        double highCutoff;
        double low;
        double high;
        NormalDistribution normal;
        for (int i = 0; i < temp.length; i++) {
            highCutoff = Math.sqrt(prior.getCutoff().getParameterValue(row * LFM.getLoadings().getColumnDimension() + i));
            lowCutoff = -highCutoff;
            for (int j = 0; j <= i; j++) {
//                if(Cholesky[i][i] > 0) {
                    if (i != j) {
                        lowCutoff = lowCutoff - temp[j] * Cholesky[i][j];
                        highCutoff = highCutoff - temp[j] * Cholesky[i][j];
                    } else {
                        lowCutoff = lowCutoff / Cholesky[i][j];
                        highCutoff = highCutoff / Cholesky[i][j];
                    }
//                }
//                else{
//                    if (i != j) {
//                        cutoffs = cutoffs + temp[j] * Cholesky[i][j];
//                    } else {
//                        cutoffs = cutoffs / Cholesky[i][j];
//                    }
//                }
            }
//            System.out.println(cutoffs);
            normal = new NormalDistribution(mean[i], 1);
            low = normal.cdf(lowCutoff);
            high = normal.cdf(highCutoff);
//            System.out.println("low: " + low);
//            System.out.println("high: " + high);
            double proportion = low/(low + 1 - high);
            if(Random.nextDouble()<proportion){
                double quantile=Random.nextDouble() * low;
                temp[i] = normal.quantile(quantile);
            }
            else{
                double quantile=(1-high) * Random.nextDouble() + high;
                temp[i] = normal.quantile(quantile);
            }

        }
        for (int i = 0; i <mean.length ; i++) {
            for (int j = 0; j <= i; j++) {
                draws[i] += Cholesky[i][j] * temp[j];
//                System.out.println("temp: " + temp[i]);
//                System.out.println("Cholesky " + i + ", " + j +": " +Cholesky[i][j]);
            }
            if(Math.abs(draws[i])<Math.sqrt(prior.getCutoff().getParameterValue(row * LFM.getLoadings().getColumnDimension() + i))) {
                System.out.println(Math.sqrt(prior.getCutoff().getParameterValue(row * LFM.getLoadings().getColumnDimension() + i)));
                System.out.println("draws: " + draws[i]);
            }
        }


        return draws;
    }

    private void drawI(int i, ListIterator<double[][]> currentPrecision, ListIterator<double[]> currentMidMean, ListIterator<double[]> currentMean) {
        double[] draws = null;
        double[][] precision = null;
        double[][] variance;
        double[] midMean = null;
        double[] mean = null;
        double[][] cholesky = null;
        if (currentPrecision.hasNext()) {
            precision = currentPrecision.next();
        }
        else
            precision = currentPrecision.previous();

        if (currentMidMean.hasNext()) {
            midMean = currentMidMean.next();
        }
        else
            midMean = currentMidMean.previous();
        if (currentMean.hasNext()) {
            mean = currentMean.next();
        }
        else
            mean= currentMean.previous();
        getPrecision(i, precision);
        variance = (new SymmetricMatrix(precision)).inverse().toComponents();

        try {
            cholesky = new CholeskyDecomposition(variance).getL();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        getMean(i, variance, midMean, mean);


        draws = getDraws(i, mean, cholesky);
//    if(i<draws.length)
//
//    {
//        while (draws[i] < 0) {
//            draws = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);
//        }
//    }
        if (i < draws.length) {
            //if (draws[i] > 0) { TODO implement as option
            copy(i, draws);
            //}
        } else {
            copy(i, draws);
        }

//       copy(i, draws);

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
        if (!randomScan) {
            ListIterator<double[][]> currentPrecision = precisionArray.listIterator();
            ListIterator<double[]> currentMidMean = meanMidArray.listIterator();
            ListIterator<double[]> currentMean = meanArray.listIterator();
            for (int i = 0; i < size; i++) {
                drawI(i, currentPrecision, currentMidMean, currentMean);
            }
            ((Parameter) loadings).fireParameterChangedEvent();
        } else {
            int i = MathUtils.nextInt(LFM.getLoadings().getRowDimension());
            ListIterator<double[][]> currentPrecision;
            ListIterator<double[]> currentMidMean;
            ListIterator<double[]> currentMean;
            if (i < LFM.getFactorDimension()) {
                currentPrecision = precisionArray.listIterator(LFM.getFactorDimension() - i - 1);
                currentMidMean = meanMidArray.listIterator(LFM.getFactorDimension() - i - 1);
                currentMean = meanArray.listIterator(LFM.getFactorDimension() - i - 1);
            } else {
                currentPrecision = precisionArray.listIterator();
                currentMidMean = meanMidArray.listIterator();
                currentMean = meanArray.listIterator();
            }
            drawI(i, currentPrecision, currentMidMean, currentMean);
            LFM.getLoadings().fireParameterChangedEvent(i, null);
//            LFM.getLoadings().fireParameterChangedEvent();
        }
        return 0;
    }

    public void setPathParameter(double beta){
        pathParameter=beta;
    }
}

