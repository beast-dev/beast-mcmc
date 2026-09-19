package dr.evomodel.coalescent;

import dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator;
import dr.inference.model.*;
import dr.math.MathUtils;
import no.uib.cipr.matrix.*;
import java.util.logging.Logger;

public class GMRFGaussianApproximation implements ModelListener, VariableListener {

    private final GMRFMultilocusSkyrideLikelihood likelihood;
    private final double tolerance;
    private final int maxNewtonIterations;
    private final int numImportanceSamples;

    // mode for approximation
    private double[] gammaHat;
    // coalescent likelihood can be approximated by a Gaussian density
    // with variable equal to the log effective population size vector
    // We denote the mean and variance of this Gaussian by \mu_g and \Sigma_g

    // Diagonal of \Sigma_g^{-1}
    private double[] sigmaGInverseDiag;
    // \Sigma_g^{-1} \mu_g
    private double[] sigmaGInverseMuG;
    private SymmTridiagMatrix scaledQ;
    private double[][] designMatrixZ;

    private boolean dirty = true;
    // to track non-convergence of Newton-Raphson method
    private boolean lastRecomputeFailed = false;

    // importance samples
    private double[][] isSamples;
    private double[] isNormalizedWeights;
    private double isESS;
    private double[] isMeanGamma;
    private boolean importanceSamplesDirty = true;

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel.coalescent.GMRFGaussianApproximation");

    public GMRFGaussianApproximation(GMRFMultilocusSkyrideLikelihood likelihood,
                                     double tolerance, int maxNewtonIterations,
                                     int numImportanceSamples) {

        this.likelihood = likelihood;
        this.tolerance = tolerance;
        this.maxNewtonIterations = maxNewtonIterations;
        this.numImportanceSamples = numImportanceSamples;

        likelihood.addModelListener(this);
        likelihood.getPrecisionParameter().addParameterListener(this);
        for(Parameter b : likelihood.getBetaListParameter()){
            b.addParameterListener(this);
        }
        likelihood.getLambdaParameter().addParameterListener(this);

        // In case we have unobserved covariate values
        for (MatrixParameter cov : likelihood.getCovariates()){
           cov.addParameterListener(this);
        }

        Parameter popSizeParameter = likelihood.getPopSizeParameter();
        this.gammaHat = new double[popSizeParameter.getDimension()];
        for(int i = 0; i < popSizeParameter.getDimension(); i++){
            gammaHat[i] = popSizeParameter.getParameterValue(i);
        }
    }

    public void modelChangedEvent(Model model, Object object, int index){
        dirty = true;
    }

    public void modelRestored(Model model){
        dirty = true;
    }

    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        dirty = true;
    }

    // Recompute gammaHat and other quantities in case
    // anything has changed since the last call
    public void ensureValuesUpToDate(){
        if(!dirty){
            return;
        }

        double tau = likelihood.getPrecisionParameter().getParameterValue(0);
        double lambda = likelihood.getLambdaParameter().getParameterValue(0);

        scaledQ = likelihood.getScaledWeightMatrix(tau, lambda);
        double[] c = likelihood.getNumCoalEvents();
        double[] suffStat = likelihood.getSufficientStatistics();
        double[] zBeta = likelihood.getZBetaVect();
        designMatrixZ = likelihood.getDesignMat();

        DenseVector gammaHatVector = new DenseVector(gammaHat);
        DenseVector mode = GMRFMultilocusSkyrideBlockUpdateOperator.newNewtonRaphson(
                c, suffStat, gammaHatVector, scaledQ.copy(),
                maxNewtonIterations, tolerance, new DenseVector(zBeta));

        if (mode == null){
            lastRecomputeFailed = true;
            LOGGER.warning("GMRFGaussianApproximation: Newton-Raphson did not converge." +
                    "We are retaining the previous gammaHat. Downstream statistics may be stale until" +
                    "we resolve this on a later iteration. Consider increasing maxNewtonIterations.");
            dirty = false;
            return;
        }
        lastRecomputeFailed = false;

        gammaHat = new double[mode.size()];
        for (int i = 0; i < gammaHat.length; i++){
            gammaHat[i] = mode.get(i);
        }

        sigmaGInverseDiag = new double[gammaHat.length];
        sigmaGInverseMuG = new double[gammaHat.length];
        for (int k = 0; k < gammaHat.length; k++){
            double expTerm = suffStat[k] * Math.exp(-gammaHat[k]);
            sigmaGInverseDiag[k] = expTerm;
            sigmaGInverseMuG[k] = -c[k] + expTerm + expTerm*gammaHat[k];
        }
        dirty = false;
    }

    // returns upper triangular band of scaledQ + sigmaGInverseDiag
    // for reuse
    private UpperTriangBandMatrix factorSystem(){
        int n = gammaHat.length;
        UpperSPDBandMatrix aMat = new UpperSPDBandMatrix(scaledQ, 1);
        for (int i = 0; i < n; i++){
            aMat.set(i, i, aMat.get(i,i)+sigmaGInverseDiag[i]);
        }
        BandCholesky cholesky = new BandCholesky(n, 1, true);
        cholesky.factor(aMat.copy());
        return cholesky.getU();
    }

    // solves (scaledQ + sigmaGInverseDiag)*x = b for each column b of bMat
    public double[][] solveSystem(double[][] bMat){
        ensureValuesUpToDate();

        int n = gammaHat.length;
        int p = bMat[0].length;
        UpperTriangBandMatrix uMat = factorSystem();
        double[][] xMat = new double[n][p];
        DenseVector tempValue = new DenseVector(n);
        DenseVector bCol = new DenseVector(n);
        DenseVector xCol = new DenseVector(n);
        for (int j = 0; j < p; j++){
            for (int i = 0; i < n; i++){
                bCol.set(i, bMat[i][j]);
            }
            uMat.transSolve(bCol, tempValue);
            uMat.solve(tempValue, xCol);
            for (int i =0; i < n; i++){
                xMat[i][j] = xCol.get(i);
            }
        }
        return xMat;
    }

    // same as above method, but with vector b rather than matrix bMat
    public double[] solveSystemVector(double[] b){
        double[][] bMat = new double[b.length][1];
        for(int i = 0; i < b.length; i++){
            bMat[i][0] = b[i];
        }
        double[][] xMat = solveSystem(bMat);
        double[] x = new double[b.length];
        for(int i = 0; i < b.length; i++){
            x[i] = xMat[i][0];
        }
        return x;
    }

    // computes diagonal of (scaledQ + sigmaGInverseDiag)^{-1}
    public double[] getGammaTildeVarianceDiag(){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        UpperTriangBandMatrix uMat = factorSystem();

        double[] varDiag = new double[n];
        DenseVector e = new DenseVector(n);
        DenseVector tempValue = new DenseVector(n);
        DenseVector x = new DenseVector(n);
        for (int i = 0; i < n; i++){
            e.zero();
            e.set(i, 1.0);
            uMat.transSolve(e, tempValue);
            uMat.solve(tempValue, x);
            varDiag[i] = x.get(i);
        }
        return varDiag;
    }

    // computes (scaledQ + sigmaGInverseDiag)^{-1} (scaledQ*Z*beta + sigmaGInverseMuG)
    // scaledQ is tau*Q
    public double[] getGammaTildeMean(){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        double[] zBeta = likelihood.getZBetaVect();
        DenseVector zBetaVector = new DenseVector(zBeta);
        DenseVector tauQZBeta = new DenseVector(n);
        scaledQ.mult(zBetaVector, tauQZBeta);

        double[] val = new double[n];
        for (int i = 0; i < n; i++){
            val[i] = tauQZBeta.get(i) + sigmaGInverseMuG[i];
        }
        return solveSystemVector(val);
    }

    public void ensureImportanceSamplesUpToDate(){
        ensureValuesUpToDate();
        if(!importanceSamplesDirty){
            return;
        }

        int n = gammaHat.length;
        double[] tildeMean = getGammaTildeMean();
        UpperTriangBandMatrix uMat = factorSystem();
        double[] c = likelihood.getNumCoalEvents();
        double[] suffStat = likelihood.getSufficientStatistics();

        isSamples = new double[numImportanceSamples][n];
        double[] logWeights = new double[numImportanceSamples];

        DenseVector z = new DenseVector(n);
        DenseVector zTransf = new DenseVector(n);

        for(int s = 0; s < numImportanceSamples; s++) {
            for (int i = 0; i < n; i++) {
                z.set(i, MathUtils.nextGaussian());
            }
            // make sure samples have appropriate variance
            uMat.solve(z, zTransf);

            double[] gammaStar = new double[n];
            double logExact = 0.0;
            double logApprox = 0.0;
            for (int i = 0; i < n; i++) {
                gammaStar[i] = tildeMean[i] + zTransf.get(i);
                logExact += -gammaStar[i] * c[i] - suffStat[i] * Math.exp(-gammaStar[i]);
                logApprox += -0.5 * sigmaGInverseDiag[i] * gammaStar[i] * gammaStar[i]
                        + sigmaGInverseMuG[i] * gammaStar[i];
            }
            isSamples[s] = gammaStar;
            logWeights[s] = logExact - logApprox;
        }

        // Normalize, using log-sum-exp trick for numerical stability
        double maxLogWeight = Double.NEGATIVE_INFINITY;
        for (int s = 0; s < numImportanceSamples; s++) {
            if(logWeights[s] > maxLogWeight) {
                maxLogWeight = logWeights[s];
            }
        }
        double sumExp = 0.0;
        double[] rawWeights = new double[numImportanceSamples];
        double sumSquaredWeights = 0.0;
        for (int s = 0; s < numImportanceSamples; s++) {
            isNormalizedWeights[s] = rawWeights[s]/sumExp;
            sumSquaredWeights += isNormalizedWeights[s]*isNormalizedWeights[s];
        }
        // Effective sample size (Kish, 1965)
        isESS = 1.0/sumSquaredWeights;

        isMeanGamma = new double[n];
        for(int s = 0; s < numImportanceSamples; s++) {
            double w = isNormalizedWeights[s];
            for(int i = 0; i < n; i++) {
                isMeanGamma[i] += w*isSamples[s][i];
            }
        }
        importanceSamplesDirty = false;
    }

    public double getImportanceSamplesESS(){
        ensureImportanceSamplesUpToDate();
        return isESS;
    }

    public int getNumImportanceSamples(){
        return numImportanceSamples;
    }

    // Importance-sampling based estimate of "information loss" term in
    // exact Fisher information J(\beta|\tau).
    // Computed as a weighted sum of P-dim rank-one outer products
    // \sum_{i} w_i*u_i*u'_i, u_i = Z'*scaledQ*(\gamma^{*(i)}-isMeanGamma)
    // this is equal to \tau^{2}*Z'QVar_{IS}[\gamma]QZ
    public double[][] getISLossMatrix(){
        ensureImportanceSamplesUpToDate();
        int n = gammaHat.length;
        int p = designMatrixZ[0].length;

        double[][] loss = new double[p][p];
        DenseVector diffVec = new DenseVector(n);
        DenseVector qDiff = new DenseVector(n);
        double[] qDiffArray = new double[n];

        for(int s = 0; s < numImportanceSamples; s++) {
            for (int i = 0; i < n; i++) {
                diffVec.set(i, isSamples[s][i]-isMeanGamma[i]);
            }
            scaledQ.mult(diffVec, qDiff);
            for (int i = 0; i < n; i++) {
                qDiffArray[i] = qDiff.get(i);
            }
            double[] u = GMRFDenseMatrixUtils.transposeMultiplyVector(designMatrixZ, qDiffArray);
            double w = isNormalizedWeights[s];
            for(int a = 0; a < p; a++){
                for (int b = 0; b < p; b++){
                    loss[a][b] += w*u[a]*u[b];
                }
            }
        }
        return loss;
    }

    // returns Q without it being scaled by tau
    public SymmTridiagMatrix getRawQ(){
        double lambda = likelihood.getLambdaParameter().getParameterValue(0);
        return likelihood.getScaledWeightMatrix(1.0, lambda);
    }

    public double[] getGammaHat(){
        ensureValuesUpToDate();
        return gammaHat;
    }

    public double[] getSigmaGInverseDiag(){
        ensureValuesUpToDate();
        return sigmaGInverseDiag;
    }

    public double[] getSigmaGInverseMuG(){
        ensureValuesUpToDate();
        return sigmaGInverseMuG;
    }

    public double[] getSuffStat(){
        return likelihood.getSufficientStatistics();
    }

    public double getPrecision(){
        return likelihood.getPrecisionParameter().getParameterValue(0);
    }

    public double[][] getDesignMatrixZ(){
        ensureValuesUpToDate();
        return designMatrixZ;
    }

    public SymmTridiagMatrix getScaledQ(){
        ensureValuesUpToDate();
        return scaledQ;
    }

    public GMRFMultilocusSkyrideLikelihood getLikelihood(){
        return likelihood;
    }

    public boolean isLastRecomputeFailed() {
        return lastRecomputeFailed;
    }
}
