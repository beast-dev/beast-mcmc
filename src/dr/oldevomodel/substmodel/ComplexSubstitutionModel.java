/*
 * ComplexSubstitutionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.oldevomodel.substmodel;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Property;
import dr.evolution.datatype.DataType;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.RobustSingularValueDecomposition;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * <b>A general irreversible class for any
 * data type; allows complex eigenstructures.</b>
 *
 * @author Marc Suchard
 */

public class ComplexSubstitutionModel extends AbstractSubstitutionModel implements Likelihood, Citable {

    public ComplexSubstitutionModel(String name, DataType dataType,
                                    FrequencyModel rootFreqModel, Parameter parameter) {

        super(name, dataType, rootFreqModel);
        this.infinitesimalRates = parameter;

        rateCount = stateCount * (stateCount - 1);

        if (parameter != null) {
            if (rateCount != infinitesimalRates.getDimension()) {
                throw new RuntimeException("Dimension of '" + infinitesimalRates.getId() + "' ("
                        + infinitesimalRates.getDimension() + ") must equal " + rateCount);
            }
            addVariable(infinitesimalRates);
        }


        stationaryDistribution = new double[stateCount];
        storedStationaryDistribution = new double[stateCount];

    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == freqModel)
            return; // freqModel only affects the likelihood calculation at the tree root
        super.handleModelChangedEvent(model, object, index);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if (!updateMatrix) {
        updateMatrix = true;
//            fireModelChanged();
//        }
    }

    protected void restoreState() {

        // To restore all this stuff just swap the pointers...

        double[] tmp3 = storedEvalImag;
        storedEvalImag = EvalImag;
        EvalImag = tmp3;

        tmp3 = storedStationaryDistribution;
        storedStationaryDistribution = stationaryDistribution;
        stationaryDistribution = tmp3;

//        normalization = storedNormalization;

        // Inherited
        updateMatrix = storedUpdateMatrix;
        wellConditioned = storedWellConditioned;

        double[] tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;

    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;

//        if(updateMatrix)
//            System.err.println("Storing updatable state!");

        storedWellConditioned = wellConditioned;

        System.arraycopy(stationaryDistribution, 0, storedStationaryDistribution, 0, stateCount);
        System.arraycopy(EvalImag, 0, storedEvalImag, 0, stateCount);
//        storedNormalization = normalization;

        // Inherited
        System.arraycopy(Eval, 0, storedEval, 0, stateCount);
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, stateCount);
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, stateCount);
        }
    }

    public void getTransitionProbabilities(double distance, double[] matrix) {

        double temp;

        int i, j, k;

        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        if (!wellConditioned) {
            Arrays.fill(matrix, 0.0);
            return;
        }


// Eigenvalues and eigenvectors of a real matrix A.
//
// If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
// and the eigenvector matrix V is orthogonal. I.e. A = V D V^t and V V^t equals
// the identity matrix.
//
// If A is not symmetric, then the eigenvalue matrix D is block diagonal with
// the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
// lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns
// of V represent the eigenvectors in the sense that A*V = V*D. The matrix
// V may be badly conditioned, or even singular, so the validity of the
// equation A = V D V^{-1} depends on the conditioning of V.

        double[][] iexp = popiexp();

        for (i = 0; i < stateCount; i++) {

            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(distance * Eval[i]);
                for (j = 0; j < stateCount; j++) {
                    iexp[i][j] = Ievc[i][j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (j = 0; j < stateCount; j++) {
                    iexp[i][j] = expatcosbt * Ievc[i][j] + expatsinbt * Ievc[i2][j];
                    iexp[i2][j] = expatcosbt * Ievc[i2][j] - expatsinbt * Ievc[i][j];
                }
                i++; // processed two conjugate rows
            }
        }

        int u = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (k = 0; k < stateCount; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }
                if (temp < 0.0)
                    matrix[u] = minProb;
                else
                    matrix[u] = temp;
                u++;
            }
        }
        pushiexp(iexp);
    }

    public double[] getStationaryDistribution() {
        return stationaryDistribution;
    }

    protected void computeStationaryDistribution() {
        stationaryDistribution = freqModel.getFrequencies();
    }


    protected double[] getRates() {
        return infinitesimalRates.getParameterValues();
    }

    protected double[] getPi() {
        return freqModel.getFrequencies();
    }

    public void setupMatrix() {

        if (!eigenInitialised) {
            initialiseEigen();
            storedEvalImag = new double[stateCount];
        }

        int i = 0;

        storeIntoAmat();


        makeValid(amat, stateCount);

        // compute eigenvalues and eigenvectors
//        EigenvalueDecomposition eigenDecomp = new EigenvalueDecomposition(new DenseDoubleMatrix2D(amat));

        RobustEigenDecomposition eigenDecomp;
        try {
            eigenDecomp = new RobustEigenDecomposition(new DenseDoubleMatrix2D(amat), maxIterations);
        } catch (ArithmeticException ae) {
            System.err.println(ae.getMessage());
            wellConditioned = false;
            System.err.println("amat = \n" + new Matrix(amat));
            return;
        }

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
        DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();
        DoubleMatrix2D eigenVInv;

        // A better (?) approach to checking diagonalizability comes from:
        //
        // J. Gentle (2007) Matrix Algebra
        //
        // Diagonalizbility Theorem: A matrix A is (complex) diagonalizable iff all distinct eigenvalues \lambda_l
        // with algebraic multiplicity m_l are semi-simple, i.e.
        //
        //          rank(A - \lambda_l I) = n - m_l
        //
        // Equivalently (?), eigenV must be non-singular.
        //
        // SVD is needed to numerically approximate the rank of a matrix, so we can check Algrebra.rank()
        // or Algebra.cond() with almost equal amounts of work.  I don't know which is more reliable. -- MAS

        if (checkConditioning) {
            RobustSingularValueDecomposition svd;
            try {
                svd = new RobustSingularValueDecomposition(eigenV, maxIterations);
            } catch (ArithmeticException ae) {
                System.err.println(ae.getMessage());
                wellConditioned = false;
                return;
            }
            if (svd.cond() > maxConditionNumber) {
                wellConditioned = false;
                return;
            }
        }

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            wellConditioned = false;
            return;
        }

        Ievc = eigenVInv.toArray();
        Evec = eigenV.toArray();
        Eval = eigenVReal.toArray();
        EvalImag = eigenVImag.toArray();

        // Check for valid decomposition
        for (i = 0; i < stateCount; i++) {
            if (Double.isNaN(Eval[i]) || Double.isNaN(EvalImag[i]) ||
                    Double.isInfinite(Eval[i]) || Double.isInfinite(EvalImag[i])) {
                wellConditioned = false;
                return;
            } else if (Math.abs(Eval[i]) < 1e-10) {
                Eval[i] = 0.0;
            }
        }

        updateMatrix = false;
        wellConditioned = true;
        // compute normalization and rescale eigenvalues

        computeStationaryDistribution();

        if (doNormalization) {
            double subst = 0.0;

            for (i = 0; i < stateCount; i++)
                subst += -amat[i][i] * stationaryDistribution[i];

//        normalization = subst;

            for (i = 0; i < stateCount; i++) {
                Eval[i] /= subst;
                EvalImag[i] /= subst;
            }
        }
    }

    //store the infinitesimal rates in the vector to a matrix called amat
    //store the infinitesimal rates in the vector to a matrix called amat
    public void storeIntoAmat(){
        double[] rates = getRates();
        double[] pi = getPi();
        int i, j, k = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = i+1; j < stateCount; j++) {
                amat[i][j] = rates[k++] * pi[j];
            }
        }
        // Copy lower triangle in column-order form (transposed)
        for (j = 0; j< stateCount; j++) {
            for (i = j+1; i < stateCount; i++) {
                amat[i][j] = rates[k++] * pi[j];
            }
        }
    }

    private void printDebugSetupMatrix() {
        System.out.println("Normalized infinitesimal rate matrix:");
        System.out.println(new Matrix(amat));
        System.out.println(new Matrix(amat).toStringOctave());
//        System.out.println("Normalization = " + normalization);
        System.out.println("Values in setupMatrix():");
//		System.out.println(eigenV);
//		System.out.println(eigenVInv);
//		System.out.println(eigenVReal);
    }

    protected void checkComplexSolutions() {
        boolean complex = false;
        for (int i = 0; i < stateCount && !complex; i++) {
            if (EvalImag[i] != 0)
                complex = true;
        }
        isComplex = complex;
    }

    public boolean getIsComplex() {
        return isComplex;
    }

    protected void frequenciesChanged() {
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates() {
    }

    protected Parameter infinitesimalRates;

//    public LogColumn[] getColumns() {
//
//        LogColumn[] columnList = new LogColumn[stateCount * stateCount];
//        int index = 0;
//        for (int i = 0; i < stateCount; i++) {
//            for (int j = 0; j < stateCount; j++)
//                columnList[index++] = new MatrixEntryColumn(getId(), i, j, amat);
//        }
//        return columnList;
//    }

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }


    public static void main(String[] arg) {

//        Parameter rates = new Parameter.Default(new double[]{5.0, 1.0, 1.0, 0.1, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
//		Parameter rates = new Parameter.Default(new double[] {5.0, 1.0, 1.0, 1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
//		Parameter rates = new Parameter.Default(new double[] {1.0, 1.0});


        Parameter rates = new Parameter.Default(159600, 1.0);


        DataType dataType = new DataType() {

            public String getDescription() {
                return null;
            }

            public int getType() {
                return 0;
            }

            @Override
            public char[] getValidChars() {
                return null;
            }

            public int getStateCount() {
                return 400;
            }
        };

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(400, 1.0 / 400.0));

        ComplexSubstitutionModel substModel = new ComplexSubstitutionModel("test",
//				TwoStates.INSTANCE,
                dataType,
                freqModel,
                rates);

        long start = System.currentTimeMillis();
        double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
        double time = 1.0;
        substModel.getTransitionProbabilities(time, finiteTimeProbs);
        long end = System.currentTimeMillis();
        System.out.println("Time: " + (end - start));
//        System.out.println("Results:");
//        System.out.println(new Vector(finiteTimeProbs));

//		System.out.println("COLT value:");
//		 This should work, matches 'octave' results
//		DoubleMatrix2D result = alegbra.mult(substModel.eigenV, alegbra.mult(blockDiagonalExponential(1.0, substModel.eigenD), substModel.eigenVInv));
//
//		System.out.println(result);

    }

    public void setMaxIterations(int max) {
        maxIterations = max;
    }

    public void setMaxConditionNumber(double max) {
        maxConditionNumber = max;
    }

    public void setCheckConditioning(boolean check) {
        checkConditioning = check;
    }

    private static DoubleMatrix2D blockDiagonalExponential(double distance, DoubleMatrix2D mat) {
        for (int i = 0; i < mat.rows(); i++) {
            if ((i + 1) < mat.rows() && mat.getQuick(i, i + 1) != 0) {
                double a = mat.getQuick(i, i);
                double b = mat.getQuick(i, i + 1);
                double expat = Math.exp(distance * a);
                double cosbt = Math.cos(distance * b);
                double sinbt = Math.sin(distance * b);
                mat.setQuick(i, i, expat * cosbt);
                mat.setQuick(i + 1, i + 1, expat * cosbt);
                mat.setQuick(i, i + 1, expat * sinbt);
                mat.setQuick(i + 1, i, -expat * sinbt);
                i++; // processed two entries in loop
            } else
                mat.setQuick(i, i, Math.exp(distance * mat.getQuick(i, i))); // 1x1 block
        }
        return mat;
    }

    private boolean isComplex = false;
    private double[] stationaryDistribution = null;
    private double[] storedStationaryDistribution;

    protected boolean doNormalization = true;
//    private Double normalization;
//    private Double storedNormalization;

    protected double[] EvalImag;
    protected double[] storedEvalImag;

    protected boolean wellConditioned = true;
    private boolean storedWellConditioned;
//    private double[] illConditionedProbabilities;

    protected static final double minProb = Property.DEFAULT.tolerance();
    //    private static final double minProb = 1E-20;
    //    private static final double minProb = Property.ZERO.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);
    EigenvalueDecomposition eigenDecomp;
    EigenvalueDecomposition storedEigenDecomp;

    private double maxConditionNumber = 1000;

    private int maxIterations = 1000;

    private boolean checkConditioning = true;

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability,this))
            return 0;
        return Double.NEGATIVE_INFINITY;
    }

    /**
     * Needs to be evaluated before the corresponding data likelihood.
     * @return
     */
    public boolean evaluateEarly() {
        return true;
    }

    public String prettyName() {
        return Abstract.getPrettyName(this);
    }

    public void setNormalization(boolean doNormalization) {
        this.doNormalization = doNormalization;
    }

    public void makeDirty() {

    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
    }

    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    private double[] probability = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Complex-diagonalizable, irreversible substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.EDWARDS_2011_ANCIENT);
    }

}