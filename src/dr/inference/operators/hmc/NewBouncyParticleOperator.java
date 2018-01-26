/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;
import java.util.Arrays;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.preorder.TipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.TipGradientViaFullConditionalDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.KroneckerOperation;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class NewBouncyParticleOperator extends SimpleMCMCOperator implements GibbsOperator {

    private double t = 0.05; //todo randomize the length a little bit.
//    private double[] v;
//    private double[] location;
//    private double[] phi_w;
//    private double[] mu;

    final NormalDistribution drawDistribution;

//    final GradientWrtParameterProvider gradientProvider;
//     double[][] precisionMatrix;
//     double[][] pMatrix;
//     double[][] sigma0; //np*np

    // Something that returns the conditional distribution
    protected Parameter parameter;

    private final TreeDataLikelihood treeDataLikelihood;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    private final TreeTrait gradientProvider; //already the gradient, not really a provider.
    private final TreeTrait densityProvider;
    private  Tree tree;

    //gradient/densityprovider is for each taxa? so precisionMatrix not necessary?( densityprovider is enough)  divide mu and p matrix from the compact form?
    // sigma0 is matrix np by np?



    public NewBouncyParticleOperator(CoercionMode mode, double weight,
                                     //TreeDataLikelihood treeDataLikelihood,
                                     ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                     String traitName,Parameter parameter, double drawVariance) {

        setWeight(weight);

        this.treeDataLikelihood  = likelihoodDelegate.getCallbackLikelihood();
        this.likelihoodDelegate = likelihoodDelegate;
        this.parameter = parameter;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));
        this.tree = treeDataLikelihood.getTree();

        String gradientName = TipGradientViaFullConditionalDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(gradientName) == null) {
            likelihoodDelegate.addFullConditionalGradientTrait(traitName);
        }
        gradientProvider = treeDataLikelihood.getTreeTrait(gradientName);
        assert (gradientProvider != null);

        String fcdName = TipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addFullConditionalDensityTrait(traitName);
        }
        densityProvider = treeDataLikelihood.getTreeTrait(fcdName);

        assert (densityProvider != null);


//        phi_w = getPhiw();
//        mu = getMU();

    }

    private double[][] getTreeVariance() {


        double priorSampleSize = likelihoodDelegate.getRootProcessDelegate().getPseudoObservations();

        return MultivariateTraitDebugUtilities.getTreeVariance(tree, treeDataLikelihood.getBranchRateModel(), 1.0,
                /*Double.POSITIVE_INFINITY*/ priorSampleSize);
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy Particle operator";
    }

    private long count = 0;

    private static final boolean DEBUG = false;

    @Override
    public double doOperation() {


        bpsOneStep();

        return 0.0;
    }

    public void setParameter(double[] position) {

        final int dim = position.length;
        for (int j = 0; j < dim; ++j) {
            parameter.setParameterValueQuietly(j, position[j]);
        }
        parameter.fireParameterChangedEvent();  // Does not seem to work with MaskedParameter
//                parameter.setParameterValueNotifyChangedAll(0, position[0]);
    }

    public void bpsOneStep() {


        double[] location = getInitialPosition();
        //System.err.println("getInitialPosition +  " + Arrays.toString(location));

        double[][] pMatrix = likelihoodDelegate.getPrecisionParameter().getParameterAsMatrix();
        double[][] sigma0 = getTreeVariance();
//        double[][] precisionMatrix = KroneckerOperation.product(pMatrix, sigma0);

        Matrix treeV = new Matrix(sigma0);
        double[][] treePrecision = (treeV.inverse()).toComponents();

//        DenseMatrix64F

        double[][] precisionMatrix = KroneckerOperation.product(treePrecision, pMatrix);





        double[] v = drawV(drawDistribution, location.length);

        t = 0.05;


        while (t > 0) {
            double[] phi_w = getPhiw();
            double[] mu = getMU();
            double[] w = addArray(location, mu, true);
            double[] phi_v = matrixMultiplier(precisionMatrix, v);
            double w_phi_w = getDotProduct(w, phi_w);//todo multiple precision matrix should be a class.
            double v_phi_w = getDotProduct(v, phi_w);//todo use a construct to store all of the temporary values.
            double v_phi_v = getDotProduct(v, phi_v);

            double tMin = Math.max(0.0, - v_phi_w/v_phi_v);

            double U_min = energyProvider(v, phi_w, phi_v, w, tMin);


            if( Double.isNaN(v_phi_w)){

                System.exit(-99);
            }
            double bounceTime = getBounceTime(v_phi_v, v_phi_w, U_min, w_phi_w);

            NewBouncyParticleOperator.TravelTime time_to_bdry = getTimeToBoundary(location, v);

            bpsUpdate(mu, v, phi_w, location,time_to_bdry, bounceTime, t, precisionMatrix);

            //System.err.println("location is (inside) +  " + Arrays.toString(location));

        }

        //System.err.println("location is (outside) +  " + Arrays.toString(location));
    }
//
//    public double[][] getVectorizedPrecisionMatrix(){
//
//        int nTaxa = tree.getExternalNodeCount();
//        int dimTaxa = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
//
//        double[][] vectorizedPrecisionMatrix = new double[nTaxa][dimTaxa];
//
//        int offsetOutput = 0;
//
//        for (int taxon = 0; taxon < nTaxa; ++taxon) { //todo to finish. from online.
//
//            System.arraycopy(precisionMatrix, 0,vectorizedPrecisionMatrix,offsetOutput,dimTaxa);
//            offsetOutput += dimTaxa;
//
//        }
//
//        return vectorizedPrecisionMatrix;
//    }

    public double[] getMU() { //todo combine getMU and getphiw and getVectorizedPrecisionMatrix

        int nTaxa = tree.getExternalNodeCount();
        int dimTaxa = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();

        double[] mu = new double[nTaxa * dimTaxa];

        int offsetOutput = 0;

        for (int taxon = 0; taxon < nTaxa; ++taxon) {

            double[] tmp = (double[]) densityProvider.getTrait(tree, tree.getExternalNode(0));

            // tmp: [0, ..., dim - 1] == mean
            // tmp: [dim, ..., dim * (dim + 1) -1] == precision or tmp[dim] * precisionParameter

            System.arraycopy(tmp, 0, mu, offsetOutput, dimTaxa);
            offsetOutput += dimTaxa;
        }

        return mu;

    }

    public double[] getPhiw() { //todo check + or - ?

        int offsetOutput = 0;
        int nTaxa = tree.getExternalNodeCount();
        int dimTaxa = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        double[] gradient = new double[nTaxa * dimTaxa];
        //todo get rid of the code duplication to get gradient
        for (int taxon = 0; taxon < nTaxa; ++taxon) {
            double[] taxonGradient = (double[]) gradientProvider.getTrait(tree, tree.getExternalNode(taxon));
            System.arraycopy(taxonGradient, 0, gradient, offsetOutput, taxonGradient.length);
            offsetOutput += taxonGradient.length;
        }
        return gradient;
    }

    public void bpsUpdate(double[] mu, double[] v, double[] phi_w, double[] location, NewBouncyParticleOperator.TravelTime traveltime, double bouncetime, double timeRemain,double[][] precisionMatrix) {

        if (timeRemain < Math.min(traveltime.minTime, bouncetime)) { //no event

            location = addArray(location, getConstantVector(v, timeRemain), false);
            t = 0.0;


        } else if (traveltime.minTime <= bouncetime) { //against the boundary

            location = addArray(location, getConstantVector(v, traveltime.minTime), false);
            location[traveltime.minIndex] = 0.0;
            v[traveltime.minIndex] = getConstantVector(v, -1.0)[traveltime.minIndex];
            t -=  traveltime.minTime;
            phi_w = matrixMultiplier(precisionMatrix, addArray(location, mu, true)); //todo to deal with the duplicated code for updating phi_w


        } else { //against the gradient

            location = addArray(location, getConstantVector(v, bouncetime), false);
            phi_w = matrixMultiplier(precisionMatrix,addArray(location, mu, true));
            v = bounceAgainst(v, phi_w);
            t -= bouncetime;

        }
        //System.err.println("location is (inside bpsupdate) +  " + Arrays.toString(location));
        setParameter(location);
    }

    public double[] bounceAgainst(double[] v, double[] minusGrad) {

        double[] finalV;
        double[] verticalV;

        verticalV = getConstantVector(minusGrad, getDotProduct(v, minusGrad)/getDotProduct(minusGrad, minusGrad));
        finalV = addArray(v, getConstantVector(verticalV, -2.0), false);

        return finalV;

    }

    public double getBounceTime(double v_phi_v, double v_phi_w, double U_min, double w_phi_w) {

        double a = v_phi_v;
        double b = 2 * v_phi_w;
        double c = 2.0 * Math.log(1 - MathUtils.nextDouble()) - U_min + w_phi_w; //TODO CHCEK IF ROOT EXISTS

        return (- b + Math.sqrt(b * b - 4 * a * c))/2/a;

    }

    public double energyProvider(double[] v, double[] phi_w, double[] phi_v, double[] w, double t){

        return getDotProduct(addArray(w, getConstantVector(v, t), false),
                addArray(phi_w, getConstantVector(phi_v, t), false));

    }

    static double[] drawV(final NormalDistribution distribution, final int dim) {

        double[] v = new double[dim];
        for (int i = 0; i < dim; i++) {
            v[i] = (Double) distribution.nextRandom();
        }
        return v;
    }

    public double[] matrixMultiplier(double[][] A, double[] B){

        final int dim = B.length;
        double[] mResult = new double[dim];

        for (int i = 0; i < dim; i ++){
            for (int j = 0; j< dim; j ++) {
                mResult[i] += A[i][j] * B[j];
            }

        }
        return mResult;
    }

    public NewBouncyParticleOperator.TravelTime getTimeToBoundary(double[] location, double[] v) {

        double[] travelTime = elementWiseVectorDivisionAbsValue(location, v);

        int index = 0;

        double minTime = Double.MAX_VALUE;

        for (int i = 0; i < travelTime.length; i ++) {
            if (travelTime[i] != 0.0 && location[i] * v[i] < 0) {
                if (travelTime[i] < minTime) {
                    index = i;
                    minTime = travelTime[i];
                }
            }
        }

        return new NewBouncyParticleOperator.TravelTime(travelTime, minTime, index);

    }

    public double[] elementWiseVectorDivisionAbsValue(double[] location, double[] y) {

        final int dim = location.length;
        double[] z = new double[dim];
        for (int i = 0; i < dim; i++) {
            z[i] = Math.abs(location[i] / y[i]); //todo make sure b[i] != 0
        }
        return z;
    }

    public static double[] addArray(double[] a, double[] b, boolean subtract) {

        assert (a.length == b.length);
        final int dim = a.length;


        double result[] = new double[dim];
        for (int i = 0; i < dim; i++) {
            if(!subtract) {
                result[i] = a[i] + b[i];
            } else {
                result[i] = a[i] - b[i];
            }

        }

        return result;
    }

    public static double getDotProduct(double[] a, double[] b) { //todo get rid of code duplication

        assert (a.length == b.length);
        final int dim = a.length;

        double total = 0.0;
        for (int i = 0; i < dim; i++) {

            total += a[i]*b[i];
        }
        return total;
    }

    public static double[] getConstantVector(double[] a, double c){

        double[] cx = new double[a.length];
        for (int i = 0; i < a.length; i ++){
            cx[i] = a[i] * c;
        }

        return cx;
    }

    public class TravelTime {

        double[] traveltime;
        double minTime;
        int minIndex;

        public TravelTime (double[] traveltime, double minTime, int minIndex){
            this.traveltime = traveltime;
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

    public double[] getInitialPosition() {
        return parameter.getParameterValues();
    }

    public double getT() {
        //todo decompose the sigma0

        return 1.0;
    }

    // Some draft functions by MAS

    public NewBouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                     PrecisionMatrixVectorProductProvider multiplicationProvider,
                                     double weight) {

        this.newGradientProvider = gradientProvider;
        this.multiplicationProvider = multiplicationProvider;
        this.likelihood = gradientProvider.getLikelihood();
        this.parameter = gradientProvider.getParameter();

        // TODO Remove
        this.drawDistribution = new NormalDistribution(0, 1);
//        this.precisionMatrix = null;
//        this.sigma0 = null;
        this.treeDataLikelihood = null;
        this.likelihoodDelegate = null;

        this.gradientProvider = null;
        this.densityProvider = null;
        this.tree = null;
        // End Remove

        setWeight(weight);
        checkParameterBounds(parameter);
    }

    private static void checkParameterBounds(Parameter parameter) {
        // TODO
    }

    private GradientWrtParameterProvider newGradientProvider;
    private PrecisionMatrixVectorProductProvider multiplicationProvider;
    private Likelihood likelihood;

}
