package dr.evomodel.coalescent;

import dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator;
import dr.inference.model.*;
import dr.math.MathUtils;
import no.uib.cipr.matrix.*;
import java.util.logging.Logger;
import org.apache.commons.math.special.Gamma;

public class GMRFGaussianApproximation implements ModelListener, VariableListener {

    private final GMRFMultilocusSkyrideLikelihood likelihood;
    private final double tolerance;
    private final int maxNewtonIterations;
    private final int numImportanceSamples;
    private final int numQuadraturePoints;
    private final double logTauLower;
    private final double logTauUpper;
    private final int numImportanceSamplesPerNode;

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
                                     int numImportanceSamples,
                                     int numQuadraturePoints,
                                     double logTauLower,
                                     double logTauUpper,
                                     int numImportanceSamplesPerNode) {

        this.likelihood = likelihood;
        this.tolerance = tolerance;
        this.maxNewtonIterations = maxNewtonIterations;
        this.numImportanceSamples = numImportanceSamples;
        this.numQuadraturePoints = numQuadraturePoints;
        this.logTauLower = logTauLower;
        this.logTauUpper = logTauUpper;
        this.numImportanceSamplesPerNode = numImportanceSamplesPerNode;

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
        importanceSamplesDirty = true;
    }

    public void modelRestored(Model model){
        dirty = true;
        importanceSamplesDirty = true;
    }

    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        dirty = true;
        importanceSamplesDirty = true;
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
        for (int s = 0; s < numImportanceSamples; s++) {
            rawWeights[s] = Math.exp(logWeights[s] - maxLogWeight);
            sumExp += rawWeights[s];
        }
        isNormalizedWeights = new double[numImportanceSamples];
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

    private UpperTriangBandMatrix factorSystemAtTau(double tau){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        double lambda = likelihood.getLambdaParameter().getParameterValue(0);
        SymmTridiagMatrix scaledQAtTau = likelihood.getScaledWeightMatrix(tau, lambda);
        UpperSPDBandMatrix aMat = new UpperSPDBandMatrix(scaledQAtTau, 1);
        for (int i = 0; i < n; i++){
            aMat.set(i, i, aMat.get(i,i)+sigmaGInverseDiag[i]);
        }
        BandCholesky cholesky = new BandCholesky(n, 1, true);
        cholesky.factor(aMat.copy());
        return cholesky.getU();
    }

    // log of det(\tau*Q + \Sigma_g^{-1}) at an arbitrary \tau
    public double getLogDetAtTau(double tau){
        UpperTriangBandMatrix uMat = factorSystemAtTau(tau);
        double acc = 0.0;
        for (int i = 0; i < uMat.numRows(); i++){
            acc += Math.log(uMat.get(i, i));
        }
        return 2.0*acc;
    }

    // Solves (\tau*Q + \Sigma_g^{-1})*x = b for each column b of bMat, at an arbitrary \tau
    public double[][] solveSystemAtTau(double tau, double[][] bMat){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        int p = bMat[0].length;
        UpperTriangBandMatrix uMat = factorSystemAtTau(tau);
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
            for (int i = 0; i < n; i++){
                xMat[i][j] = xCol.get(i);
            }
        }
        return xMat;
    }

    private double[] solveSystemAtTauVector(double tau, double[] b){
        double[][] bMat = new double[b.length][1];
        for (int i = 0; i < b.length; i++){
            bMat[i][0] = b[i];
        }
        double[][] xMat = solveSystemAtTau(tau, bMat);
        double[] x = new double[b.length];
        for (int i = 0; i < b.length; i++){
            x[i] = xMat[i][0];
        }
        return x;
    }

    // \tilde{mu}_{\gamma|g,Z,\beta,\tau} at an arbitrary \tau (generalizes getGammaTildeMean())
    public double[] getTildeMeanAtTau(double tau){
        ensureValuesUpToDate();
        double lambda = likelihood.getLambdaParameter().getParameterValue(0);
        int n = gammaHat.length;
        double[] zBeta = likelihood.getZBetaVect();
        SymmTridiagMatrix scaledQAtTau = likelihood.getScaledWeightMatrix(tau, lambda);
        DenseVector zBetaVector = new DenseVector(zBeta);
        DenseVector tauQZBeta = new DenseVector(n);
        scaledQAtTau.mult(zBetaVector, tauQZBeta);

        double[] val = new double[n];
        for (int i = 0; i < n; i++){
            val[i] = tauQZBeta.get(i) + sigmaGInverseMuG[i];
        }
        return solveSystemAtTauVector(tau, val);
    }


    // Result of drawing and reweighting numSamples \gamma draws from P_G(\gamma|g,Z,\beta,\tau) at a given \tau
    public static class WeightedSamples {
        public final double[][] samples;          // numSamples x n
        public final double[] normalizedWeights;   // numSamples, sums to 1
        public final double logMeanRawWeight;       // log[(1/S)*sum_i exp(logWeight_i)], i.e. log of the RAW (unnormalized) mean weight
        public final double ess;
        public final double[] weightedMean;         // length n

        WeightedSamples(double[][] samples, double[] normalizedWeights, double logMeanRawWeight,
                        double ess, double[] weightedMean){
            this.samples = samples;
            this.normalizedWeights = normalizedWeights;
            this.logMeanRawWeight = logMeanRawWeight;
            this.ess = ess;
            this.weightedMean = weightedMean;
        }
    }

    public WeightedSamples drawAndWeightSamplesAtTau(double tau, int numSamples){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        double[] tildeMeanAtTau = getTildeMeanAtTau(tau);
        UpperTriangBandMatrix uMat = factorSystemAtTau(tau);
        double[] c = likelihood.getNumCoalEvents();
        double[] suffStat = likelihood.getSufficientStatistics();

        double[][] samples = new double[numSamples][n];
        double[] logWeights = new double[numSamples];

        DenseVector z = new DenseVector(n);
        DenseVector sampleDeviation = new DenseVector(n);

        for (int s = 0; s < numSamples; s++){
            for (int i = 0; i < n; i++){
                z.set(i, MathUtils.nextGaussian());
            }
            uMat.solve(z, sampleDeviation);

            double[] gammaStar = new double[n];
            double logExact = 0.0;
            double logApprox = 0.0;
            for (int i = 0; i < n; i++){
                gammaStar[i] = tildeMeanAtTau[i] + sampleDeviation.get(i);
                logExact += -gammaStar[i]*c[i] - suffStat[i]*Math.exp(-gammaStar[i]);
                logApprox += -0.5*sigmaGInverseDiag[i]*gammaStar[i]*gammaStar[i]
                        + sigmaGInverseMuG[i]*gammaStar[i];
            }
            samples[s] = gammaStar;
            logWeights[s] = logExact - logApprox;
        }

        double[] normalizedWeights = GMRFDenseMatrixUtils.normalizeLogWeights(logWeights);
        double ess = GMRFDenseMatrixUtils.effectiveSampleSize(normalizedWeights);
        double logMeanRawWeight = GMRFDenseMatrixUtils.logSumExp(logWeights) - Math.log(numSamples);

        double[] weightedMean = new double[n];
        for (int s = 0; s < numSamples; s++){
            double w = normalizedWeights[s];
            for (int i = 0; i < n; i++){
                weightedMean[i] += w * samples[s][i];
            }
        }

        return new WeightedSamples(samples, normalizedWeights, logMeanRawWeight, ess, weightedMean);
    }

    private static class NodeFactorization {
        final UpperTriangBandMatrix uMat;
        final double logDetA;
        NodeFactorization(UpperTriangBandMatrix uMat, double logDetA){
            this.uMat = uMat;
            this.logDetA = logDetA;
        }
    }

    private NodeFactorization factorAndLogDetAtTau(double tau){
        UpperTriangBandMatrix uMat = factorSystemAtTau(tau);
        double acc = 0.0;
        for (int i = 0; i < uMat.numRows(); i++){
            acc += Math.log(uMat.get(i, i));
        }
        return new NodeFactorization(uMat, 2.0*acc);
    }

    private static double[] solveGivenFactorization(UpperTriangBandMatrix uMat, double[] b){
        int n = b.length;
        DenseVector bVec = new DenseVector(b);
        DenseVector tmp = new DenseVector(n);
        DenseVector x = new DenseVector(n);
        uMat.transSolve(bVec, tmp);
        uMat.solve(tmp, x);
        double[] result = new double[n];
        for (int i = 0; i < n; i++){
            result[i] = x.get(i);
        }
        return result;
    }

    private static double dot(double[] a, double[] b){
        double s = 0.0;
        for (int i = 0; i < a.length; i++){
            s += a[i]*b[i];
        }
        return s;
    }

    private static class MarginalSetup {
        final double[] zBeta;
        final double[] rawQZBeta;     // Q_raw * Z * beta
        final double zBetaRawQZBeta;  // (Z*beta)' * Q_raw * (Z*beta)
        final double[][] rawQZ;       // Q_raw * Z, (M+1) x P
        final double[][] rawZQZ;      // Z' * Q_raw * Z, P x P (tau-independent)
        final double[] zqzBetaRaw;    // Z' * Q_raw * (Z*beta), P-vector

        MarginalSetup(SymmTridiagMatrix rawQ, double[][] designMatrixZ, double[] zBeta){
            this.zBeta = zBeta;
            int n = zBeta.length;
            int p = designMatrixZ[0].length;

            DenseVector zBetaVec = new DenseVector(zBeta);
            DenseVector rawQZBetaVec = new DenseVector(n);
            rawQ.mult(zBetaVec, rawQZBetaVec);
            this.rawQZBeta = new double[n];
            for (int i = 0; i < n; i++){
                rawQZBeta[i] = rawQZBetaVec.get(i);
            }
            this.zBetaRawQZBeta = dot(zBeta, rawQZBeta);

            this.rawQZ = new double[n][p];
            for (int j = 0; j < p; j++){
                DenseVector col = new DenseVector(n);
                for (int i = 0; i < n; i++){
                    col.set(i, designMatrixZ[i][j]);
                }
                DenseVector qCol = new DenseVector(n);
                rawQ.mult(col, qCol);
                for (int i = 0; i < n; i++){
                    rawQZ[i][j] = qCol.get(i);
                }
            }
            this.rawZQZ = GMRFDenseMatrixUtils.transposeMultiply(designMatrixZ, rawQZ);
            this.zqzBetaRaw = GMRFDenseMatrixUtils.transposeMultiplyVector(designMatrixZ, rawQZBeta);
        }
    }

    public static class MarginalFisherInfoResult {
        public final double[][] jMarginal; // P x P
        public final double eTau;           // E_tau[tau] under the tau-weighting used
        public final double minESS;         // NaN for the Gaussian version; min over nodes for the exact version

        MarginalFisherInfoResult(double[][] jMarginal, double eTau, double minESS){
            this.jMarginal = jMarginal;
            this.eTau = eTau;
            this.minESS = minESS;
        }
    }

    public MarginalFisherInfoResult getMarginalFisherInformationGaussian(double tauShape, double tauRate){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        int p = designMatrixZ[0].length;
        int M = n - 1;

        SymmTridiagMatrix rawQ = getRawQ();
        MarginalSetup setup = new MarginalSetup(rawQ, designMatrixZ, likelihood.getZBetaVect());

        double[][] nodesAndWeights = GMRFDenseMatrixUtils.simpsonNodesAndWeights(numQuadraturePoints, logTauLower, logTauUpper);
        double[] points = nodesAndWeights[0];
        double[] simpsonW = nodesAndWeights[1];
        int K = points.length;

        double[] logNodeWeight = new double[K];
        double[] tauAtNode = new double[K];
        double[][][] jAtNode = new double[K][][];
        double[][] scoreAtNode = new double[K][];

        for (int k = 0; k < K; k++){
            double uK = points[k];
            double tauK = Math.exp(uK);
            tauAtNode[k] = tauK;

            NodeFactorization nf = factorAndLogDetAtTau(tauK);

            double[] vK = new double[n];
            for (int i = 0; i < n; i++){
                vK[i] = tauK*setup.rawQZBeta[i] + sigmaGInverseMuG[i];
            }
            double[] xK = solveGivenFactorization(nf.uMat, vK); // tildeMean(tauK)
            double quadFormK = dot(vK, xK);

            double logPGgK = (M/2.0)*uK - 0.5*nf.logDetA + 0.5*quadFormK - 0.5*tauK*setup.zBetaRawQZBeta;
            double logPTauK = tauShape*Math.log(tauRate) - Gamma.logGamma(tauShape)
                    + (tauShape-1.0)*uK - tauRate*tauK;
            // Jacobian for the log-tau substitution folded in as "+uK" (matches
            // the (a+M/2-1)*u + u = (a+M/2)*u pattern used in SkygridMutualInformationAnalysis).
            logNodeWeight[k] = logPTauK + logPGgK + uK + Math.log(simpsonW[k]);

            double[][] wK = GMRFDenseMatrixUtils.scaleMult(setup.rawQZ, tauK);
            double[][] vMatK = new double[n][p];
            for (int j = 0; j < p; j++){
                double[] col = new double[n];
                for (int i = 0; i < n; i++){
                    col[i] = wK[i][j];
                }
                double[] solved = solveGivenFactorization(nf.uMat, col);
                for (int i = 0; i < n; i++){
                    vMatK[i][j] = solved[i];
                }
            }
            double[][] cK = GMRFDenseMatrixUtils.scaleMult(setup.rawZQZ, tauK);
            double[][] wtVK = GMRFDenseMatrixUtils.transposeMultiply(wK, vMatK);
            jAtNode[k] = GMRFDenseMatrixUtils.subtractMat(cK, wtVK);

            // score_tau^G(tauK) = tauK * Z'Q_raw * (tildeMean(tauK) - Z*beta)
            // -- reuses xK (already computed for the node weight above), plus
            // one cheap (M+1)-dim Q multiply and a Z' reduction; no new solve.
            double[] diffK = new double[n];
            for (int i = 0; i < n; i++){
                diffK[i] = xK[i] - setup.zBeta[i];
            }
            DenseVector diffVecK = new DenseVector(diffK);
            DenseVector rawQDiffK = new DenseVector(n);
            rawQ.mult(diffVecK, rawQDiffK);
            double[] rawQDiffArrK = new double[n];
            for (int i = 0; i < n; i++){
                rawQDiffArrK[i] = rawQDiffK.get(i);
            }
            double[] zqDiffK = GMRFDenseMatrixUtils.transposeMultiplyVector(designMatrixZ, rawQDiffArrK);
            double[] scoreK = new double[p];
            for (int j = 0; j < p; j++){
                scoreK[j] = tauK*zqDiffK[j];
            }
            scoreAtNode[k] = scoreK;
        }

        double[] w = GMRFDenseMatrixUtils.normalizeLogWeights(logNodeWeight);

        double[][] eJGaussian = new double[p][p];
        double eTau = 0.0;
        double[] scoreMean = new double[p];
        for (int k = 0; k < K; k++){
            eTau += w[k]*tauAtNode[k];
            for (int a = 0; a < p; a++){
                for (int b = 0; b < p; b++){
                    eJGaussian[a][b] += w[k]*jAtNode[k][a][b];
                }
                scoreMean[a] += w[k]*scoreAtNode[k][a];
            }
        }

        double[][] scoreVar = new double[p][p];
        for (int k = 0; k < K; k++){
            double[] diff = new double[p];
            for (int a = 0; a < p; a++){
                diff[a] = scoreAtNode[k][a] - scoreMean[a];
            }
            for (int a = 0; a < p; a++){
                for (int b = 0; b < p; b++){
                    scoreVar[a][b] += w[k]*diff[a]*diff[b];
                }
            }
        }

        double[][] jMarginal = GMRFDenseMatrixUtils.subtractMat(eJGaussian, scoreVar);

        return new MarginalFisherInfoResult(jMarginal, eTau, Double.NaN);
    }

    public MarginalFisherInfoResult getMarginalFisherInformationIS(double tauShape, double tauRate){
        ensureValuesUpToDate();
        int n = gammaHat.length;
        int p = designMatrixZ[0].length;
        int M = n - 1;

        SymmTridiagMatrix rawQ = getRawQ();
        MarginalSetup setup = new MarginalSetup(rawQ, designMatrixZ, likelihood.getZBetaVect());

        double[][] nodesAndWeights = GMRFDenseMatrixUtils.simpsonNodesAndWeights(numQuadraturePoints, logTauLower, logTauUpper);
        double[] points = nodesAndWeights[0];
        double[] simpsonW = nodesAndWeights[1];
        int K = points.length;

        double[] logNodeWeight = new double[K];
        double[] tauAtNode = new double[K];
        double[][][] jExactAtNode = new double[K][][];
        double[][] scoreAtNode = new double[K][];
        double minESS = Double.POSITIVE_INFINITY;

        for (int k = 0; k < K; k++){
            double uK = points[k];
            double tauK = Math.exp(uK);
            tauAtNode[k] = tauK;

            // Gaussian-approximation piece, needed for the IS-corrected node weight
            // (logPGgK below), exactly as in the Gaussian version.
            NodeFactorization nf = factorAndLogDetAtTau(tauK);
            double[] vK = new double[n];
            for (int i = 0; i < n; i++){
                vK[i] = tauK*setup.rawQZBeta[i] + sigmaGInverseMuG[i];
            }
            double[] xK = solveGivenFactorization(nf.uMat, vK);
            double quadFormK = dot(vK, xK);
            double logPGgK = (M/2.0)*uK - 0.5*nf.logDetA + 0.5*quadFormK - 0.5*tauK*setup.zBetaRawQZBeta;

            // Importance-sampling correction at this node: draws from
            // P_G(gamma|g,Z,beta,tauK), reweighted by P(g|gamma)/P_G(g|gamma).
            WeightedSamples ws = drawAndWeightSamplesAtTau(tauK, numImportanceSamplesPerNode);
            minESS = Math.min(minESS, ws.ess);

            double logPTauK = tauShape*Math.log(tauRate) - Gamma.logGamma(tauShape)
                    + (tauShape-1.0)*uK - tauRate*tauK;
            // True (IS-corrected) unnormalized P(tauK|g,Z,beta), in log space,
            // via P(g|Z,beta,tauK) = P_G(g|Z,beta,tauK) * [raw-mean correction].
            logNodeWeight[k] = logPTauK + logPGgK + ws.logMeanRawWeight + uK + Math.log(simpsonW[k]);

            // Exact J(beta|tauK) via the same rank-one IS construction as
            // getISLossMatrix(), parameterized by this node's tauK/samples.
            double[][] lossK = new double[p][p];
            DenseVector diffVec = new DenseVector(n);
            DenseVector qDiff = new DenseVector(n);
            double[] qDiffArray = new double[n];
            for (int s = 0; s < numImportanceSamplesPerNode; s++){
                for (int i = 0; i < n; i++){
                    diffVec.set(i, ws.samples[s][i] - ws.weightedMean[i]);
                }
                // tauK * Q_raw * diff, via the raw (unscaled) Q and manual scaling
                // by tauK, avoiding rebuilding an MTJ scaled-Q object per sample.
                DenseVector rawQDiff = new DenseVector(n);
                rawQ.mult(diffVec, rawQDiff);
                for (int i = 0; i < n; i++){
                    qDiff.set(i, tauK*rawQDiff.get(i));
                    qDiffArray[i] = qDiff.get(i);
                }
                double[] u = GMRFDenseMatrixUtils.transposeMultiplyVector(designMatrixZ, qDiffArray);
                double wgt = ws.normalizedWeights[s];
                for (int a = 0; a < p; a++){
                    for (int b = 0; b < p; b++){
                        lossK[a][b] += wgt * u[a] * u[b];
                    }
                }
            }
            double[][] completeK = GMRFDenseMatrixUtils.scaleMult(setup.rawZQZ, tauK);
            jExactAtNode[k] = GMRFDenseMatrixUtils.subtractMat(completeK, lossK);

            // Exact score_tau(tauK) = tauK*Z'Q*(weightedMean - Z*beta)
            //                       = tauK*(Z'*Q_raw*weightedMean - zqzBetaRaw).
            double[] rawQMean = new double[n];
            DenseVector meanVec = new DenseVector(ws.weightedMean);
            DenseVector rawQMeanVec = new DenseVector(n);
            rawQ.mult(meanVec, rawQMeanVec);
            for (int i = 0; i < n; i++){
                rawQMean[i] = rawQMeanVec.get(i);
            }
            double[] zqMean = GMRFDenseMatrixUtils.transposeMultiplyVector(designMatrixZ, rawQMean);
            double[] scoreK = new double[p];
            for (int j = 0; j < p; j++){
                scoreK[j] = tauK*(zqMean[j] - setup.zqzBetaRaw[j]);
            }
            scoreAtNode[k] = scoreK;
        }

        double[] w = GMRFDenseMatrixUtils.normalizeLogWeights(logNodeWeight);

        double[][] eJExact = new double[p][p];
        double eTau = 0.0;
        double[] scoreMean = new double[p];
        for (int k = 0; k < K; k++){
            eTau += w[k]*tauAtNode[k];
            for (int a = 0; a < p; a++){
                for (int b = 0; b < p; b++){
                    eJExact[a][b] += w[k]*jExactAtNode[k][a][b];
                }
                scoreMean[a] += w[k]*scoreAtNode[k][a];
            }
        }

        double[][] scoreVar = new double[p][p];
        for (int k = 0; k < K; k++){
            double[] diff = new double[p];
            for (int a = 0; a < p; a++){
                diff[a] = scoreAtNode[k][a] - scoreMean[a];
            }
            for (int a = 0; a < p; a++){
                for (int b = 0; b < p; b++){
                    scoreVar[a][b] += w[k]*diff[a]*diff[b];
                }
            }
        }

        double[][] jMarginal = GMRFDenseMatrixUtils.subtractMat(eJExact, scoreVar);

        return new MarginalFisherInfoResult(jMarginal, eTau, minESS);
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
