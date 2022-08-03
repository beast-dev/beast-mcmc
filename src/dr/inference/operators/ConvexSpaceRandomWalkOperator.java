package dr.inference.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import jebl.math.Random;
import org.ejml.data.DenseMatrix64F;

import java.util.ArrayList;

public class ConvexSpaceRandomWalkOperator extends AbstractAdaptableOperator {

    private static final boolean DEBUG = true;
    private static final double corDeflator = 0.01;
    private final ArrayList<double[]> sampleList = new ArrayList<>();


    private double window;
    private final BoundedSpace space;
    private final Parameter parameter;
    private final Parameter updateIndex;
    private final boolean ADAPTIVE_COVARIANCE = true;
    private final int burnin = 50;
    private double[] mean;
    private double[] oldMean;
    private final DenseMatrix64F cov;
    private final int dim;
    private final int varDim;
    private int iterations = 0;
    private int updates = 0;
    private int every = 10;
    private final ArrayList<Integer> varInds = new ArrayList<>();
    private double[][] cholesky;

    public static final String CONVEX_RW = "convexSpaceRandomWalkOperator";
    public static final String WINDOW_SIZE = "relativeWindowSize";

    public ConvexSpaceRandomWalkOperator(Parameter parameter, BoundedSpace space,
                                         Parameter updateIndex,
                                         double window, double weight) {
        setWeight(weight);

        this.updateIndex = updateIndex;
        this.parameter = parameter;
        this.space = space;
        this.window = window;

        this.dim = parameter.getDimension();
        for (int i = 0; i < dim; i++) {
            if (updateIndex == null || updateIndex.getParameterValue(i) == 1) {
                varInds.add(i);
            }
        }
        this.varDim = varInds.size();

        this.cov = new DenseMatrix64F(varDim, varDim);
        for (int i = 0; i < varDim; i++) {
            cov.set(i, i, 1);
        }
        cholesky = CholeskyDecomposition.execute(cov.getData(), 0, varDim);

    }


    @Override
    public double doOperation() {
        iterations++;
        double[] values = parameter.getParameterValues();
        double[] varValues = new double[varDim];
        for (int i = 0; i < varDim; i++) {
            varValues[i] = values[varInds.get(i)];
        }

        if (ADAPTIVE_COVARIANCE) {


            if (iterations == burnin) {
                sampleList.add(varValues);
                mean = varValues;
                updates++;
            } else if (iterations > burnin && iterations % every == 0) {
                sampleList.add(varValues);
                updates++;
                oldMean = mean;

                for (int i = 0; i < varDim; i++) {
                    mean[i] = ((updates - 1) * oldMean[i] + varValues[i]) / updates;
                }

                for (int i = 0; i < varDim; i++) {
                    for (int j = i; j < varDim; j++) {
                        double value = (updates - 1) * (cov.get(i, j) + oldMean[i] * oldMean[j]);
                        value += varValues[i] * varValues[j] - mean[i] * mean[j];
                        value /= updates;
                        cov.set(i, j, value);
                        cov.set(j, i, value);
                    }
                }

                DenseMatrix64F cor = new DenseMatrix64F(varDim, varDim);
                for (int i = 0; i < varDim; i++) {
                    for (int j = 0; j < varDim; j++) {
                        cor.set(i, j, cov.get(i, j) / Math.sqrt(cov.get(i, i) * cov.get(j, j)));
                        if (i == j) {
                            cor.set(i, j, cor.get(i, j) + corDeflator);
                        }
                    }
                }

                if (DEBUG) {

                    DenseMatrix64F testCov = new DenseMatrix64F(varDim, varDim);
                    double[] testMean = new double[varDim];

                    for (int i = 0; i < updates; i++) {
                        double[] valuesi = sampleList.get(i);
                        for (int j = 0; j < varDim; j++) {
                            testMean[j] += valuesi[j];
                            testCov.set(j, j, testCov.get(j, j) + valuesi[j] * valuesi[j]);
                            for (int k = (j + 1); k < varDim; k++) {
                                testCov.set(j, k, testCov.get(j, k) + valuesi[j] * valuesi[k]);
                            }
                        }
                    }

//                    System.out.print("Sum squares:");
//                    System.out.print(testCov);

                    for (int i = 0; i < varDim; i++) {
                        testMean[i] /= updates;
                    }

//                    System.out.print("Mean: ");
//                    Utils.printArray(testMean);

                    for (int i = 0; i < varDim; i++) {
                        testCov.set(i, i, testCov.get(i, i) / updates - testMean[i] * testMean[i]);
                        for (int j = (i + 1); j < varDim; j++) {
                            testCov.set(i, j, testCov.get(i, j) / updates - testMean[i] * testMean[j]);
                            testCov.set(j, i, testCov.get(i, j));
                        }
                    }


                    DenseMatrix64F testCor = new DenseMatrix64F(varDim, varDim);
                    for (int i = 0; i < varDim; i++) {
                        for (int j = 0; j < varDim; j++) {
                            testCor.set(i, j, testCov.get(i, j) / Math.sqrt(testCov.get(i, i) * testCov.get(j, j)));
                            if (i == j) {
                                testCor.set(i, j, testCor.get(i, j) + corDeflator);
                            }
                        }
                    }


//                    System.out.println("Cov:");
//                    System.out.println(testCov);
//                    System.out.println(cov);
//
//                    System.out.println("Cor");
//                    System.out.print(testCor);
//                    System.out.println(cor);
//                    System.out.println();

                    System.arraycopy(testCor.data, 0, cor.data, 0, cor.data.length); //TODO: remove after fixing
                }


                cholesky = CholeskyDecomposition.execute(cor.getData(), 0, varDim);
            }
        }

//        double[] sample = (double[]) generator.nextRandom();

        double[] sample = new double[values.length];
        double[] varSample = MultivariateNormalDistribution.nextMultivariateNormalCholesky(new double[varDim], cholesky);

        double sum = 0;
        for (int i = 0; i < varDim; i++) {
            sum += varSample[i] * varSample[i];
        }

        double norm = Math.sqrt(sum);
        for (int i = 0; i < varDim; i++) {
            varSample[i] = varSample[i] / norm;
        }

        for (int i = 0; i < varDim; i++) {
            sample[varInds.get(i)] = varSample[i];
        }

        BoundedSpace.IntersectionDistances distances;
        try {
            distances = space.distancesToBoundary(values, sample, false);
        } catch (HamiltonianMonteCarloOperator.NumericInstabilityException e) {
            throw new RuntimeException("position outside of bounded space at beginning of operator move");
        }
//        double u1 = Random.nextDouble() * distances.forwardDistance;
//        for (int i = 0; i < values.length; i++) {
//            sample[i] = values[i] + (sample[i] - values[i]) * u1;
//        }


        double t = window * Random.nextDouble();
        t = RandomWalkOperator.reflectValue(t, -distances.backwardDistance, distances.forwardDistance);

        for (int i = 0; i < values.length; i++) {
            sample[i] = values[i] - sample[i] * t;
        }

        parameter.setAllParameterValuesQuietly(sample);
        parameter.fireParameterChangedEvent();

//        double tForward = t / (distances.forwardDistance * window);
//        double forwardLogDensity = uniformProductLogPdf(tForward);
//
//        double backWardDistance = distances.backwardDistance + t;
//        double tBackward = t / (backWardDistance * window);
//        double backwardLogDensity = uniformProductLogPdf(tBackward);

        return 0;
    }

//    private double uniformProductLogPdf(double t) {
//        double density = -Math.log(t);
//        return Math.log(density);
//    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (value > 0) value = 0;
        window = Math.exp(value);
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(window);
    }

    @Override
    public double getRawParameter() {
        return window;
    }

    @Override
    public String getAdaptableParameterName() {
        return WINDOW_SIZE;
    }

    @Override
    public String getOperatorName() {
        return CONVEX_RW;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }


}
