package dr.evomodel.coalescent;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

public class SkygridSummaryStatistic extends Statistic.Abstract implements Loggable {

    public enum Type{
        RETAINED_INFORMATION_RATIO,
        MUTUAL_INFORMATION,
        LAGRANGE_BOUND,
        GAMMA_HAT,
        SIGMA_G_INVERSE_DIAG,
        GAMMA_TILDE_MEAN,
        GAMMA_TILDE_VARIANCE_DIAG,
        BETA_GIVEN_GAMMA_MEAN,
        BETA_GIVEN_GAMMA_VARIANCE_DIAG,
        BETA_TILDE_MEAN,
        BETA_TILDE_VARIANCE_DIAG,
        TAU_CONDITIONAL_SHAPE,
        TAU_CONDITIONAL_RATE
    }

    private final GMRFGaussianApproximation approximation;
    private final Type type;
    private final Integer coefficientIndex;
    private final double[][] betaPriorPrecision;
    private final double[] betaPriorMean;
    private final Double tauShape;
    private final Double tauRate;
    private final Parameter gamma;

    public SkygridSummaryStatistic(String name,
                                   GMRFGaussianApproximation approximation,
                                   Type type,
                                   Integer coefficientIndex,
                                   double[][] betaPriorPrecision,
                                   double[] betaPriorMean,
                                   Double tauShape,
                                   Double tauRate,
                                   Parameter gamma){
        super(name);
        this.approximation = approximation;
        this.type = type;
        this.coefficientIndex = coefficientIndex;
        this.betaPriorPrecision = betaPriorPrecision;
        this.betaPriorMean = betaPriorMean;
        this.tauShape = tauShape;
        this.tauRate = tauRate;
        this.gamma = gamma;
    }

    public int getDimension(){
        switch (type) {
            case RETAINED_INFORMATION_RATIO:
            case MUTUAL_INFORMATION:
            case TAU_CONDITIONAL_SHAPE:
            case TAU_CONDITIONAL_RATE:
                return 1;
            case LAGRANGE_BOUND:
            case GAMMA_HAT:
            case SIGMA_G_INVERSE_DIAG:
            case GAMMA_TILDE_MEAN:
            case GAMMA_TILDE_VARIANCE_DIAG:
                return approximation.getGammaHat().length;
            case BETA_GIVEN_GAMMA_MEAN:
            case BETA_GIVEN_GAMMA_VARIANCE_DIAG:
            case BETA_TILDE_MEAN:
            case BETA_TILDE_VARIANCE_DIAG:
                return approximation.getDesignMatrixZ()[0].length;
            default:
                throw new IllegalArgumentException("Unsupported type: " + type);
        }
    }

    public double getStatisticValue(int dimension){
        approximation.ensureValuesUpToDate();

        switch (type) {
            case RETAINED_INFORMATION_RATIO:
                return computeRetainedInformationRatio();
            case MUTUAL_INFORMATION:
                return computeMutualInformation();
            case LAGRANGE_BOUND:
                return computeLagrangeBoundPerInterval()[dimension];
            case GAMMA_HAT:
                return approximation.getGammaHat()[dimension];
            case SIGMA_G_INVERSE_DIAG:
                return approximation.getSigmaGInverseDiag()[dimension];
            case GAMMA_TILDE_MEAN:
                return approximation.getGammaTildeMean()[dimension];
            case GAMMA_TILDE_VARIANCE_DIAG:
                return approximation.getGammaTildeVarianceDiag()[dimension];
            case BETA_GIVEN_GAMMA_MEAN:
                return computeBetaGivenGammaMean()[dimension];
            case BETA_GIVEN_GAMMA_VARIANCE_DIAG:
                return computeBetaGivenGammaVarianceDiag()[dimension];
            case BETA_TILDE_MEAN:
                return computeBetaTildeMean()[dimension];
            case BETA_TILDE_VARIANCE_DIAG:
                return computeBetaTildeVarianceDiag()[dimension];
            case TAU_CONDITIONAL_SHAPE:
                return computeTauConditionalShape();
            case TAU_CONDITIONAL_RATE:
                return computeTauConditionalRate();
            default:
                throw new IllegalArgumentException("Unsupported type: " + type);
        }
    }

    // Computes (1) information matrix if we had "complete information": \tau*Z'QZ
    // and (2) approximate Fisher information
    // J_G(\beta|\tau) = "complete information" - \tau^{2}*Z'Q(\tau*Q+\Sigma_g^{-1})^{-1}QZ
    private double[][][] computeInformationMatrices(){
        double[][] zMat = approximation.getDesignMatrixZ();
        SymmTridiagMatrix scaledQ = approximation.getScaledQ();
         int n = zMat.length;
         int p = zMat[0].length;

         double[][] wMat = new double[n][p];
         for (int j = 0; j < p; j++) {
             DenseVector zCol = new DenseVector(n);
             for (int i = 0; i < n; i++) {
                 zCol.set(i, zMat[i][j]);
             }
             DenseVector wCol = new DenseVector(n);
             // W = \tau*QZ
             scaledQ.mult(zCol, wCol);
             for (int i = 0; i < n; i++) {
                 wMat[i][j] = wCol.get(i);
             }
         }
         // V = (\tau*Q + \Sigma_g^{-1})^{-1} *W
         double[][] vMat = approximation.solveSystem(wMat);
         // compute Z'\tau*Q Z
         double[][] complete = GMRFDenseMatrixUtils.transposeMultiply(zMat, wMat);
         // compute \tau^{2}*Z'Q(\tau*Q+\Sigma_g^{-1})^{-1}QZ
         double[][] loss = GMRFDenseMatrixUtils.transposeMultiply(wMat, vMat);
         double[][] jMat = GMRFDenseMatrixUtils.subtractMat(complete, loss);

         return new double[][][]{complete,jMat};
    }

    // compute \tau*Z'Qv for a vector v
    private double[] tauZTransposeQTimes(double[] v){
        SymmTridiagMatrix scaledQ = approximation.getScaledQ();
        double[][] zMat = approximation.getDesignMatrixZ();
        int n = v.length;
        DenseVector vVect = new DenseVector(v);
        DenseVector qv = new DenseVector(n);
        scaledQ.mult(vVect, qv);
        double[] qvArray = new double[n];
        for (int i = 0; i < n; i++) {
            qvArray[i] = qv.get(i);
        }
        return GMRFDenseMatrixUtils.transposeMultiplyVector(zMat, qvArray);
    }

    private double[] gammaAsArray(){
        int n = gamma.getDimension();
        double[] gammaArray = new double[n];
        for (int i = 0; i < n; i++) {
            gammaArray[i] = gamma.getParameterValue(i);
        }
        return gammaArray;
    }

    // computes \frac{v'J_G(\beta | \tau)v}{v'[\tau*Z'QZ]v}
    private double computeRetainedInformationRatio(){
        double[][][] informationMatrices = computeInformationMatrices();
        double[][] complete = informationMatrices[0];
        double[][] jMat = informationMatrices[1];
        int p = complete.length;

        if(coefficientIndex != null){
            double[] v = new double[p];
            v[coefficientIndex] = 1;
            double numerator = GMRFDenseMatrixUtils.quadraticForm(jMat, v);
            double denominator = GMRFDenseMatrixUtils.quadraticForm(complete, v);
            return numerator / denominator;
        }else{
            double traceComplete = 0;
            double traceJ = 0;
            for (int i = 0; i < p; i++) {
                traceComplete += complete[i][i];
                traceJ += jMat[i][i];
            }
            return traceJ / traceComplete;
        }
    }

    // Computes approximation of mutual information of \gamma and \beta, conditional
    // on g, Z and \tau. Depends on Gaussian approximation of coalescent likelihood.
    private double computeMutualInformation(){
        double[][][] informationMatrices = computeInformationMatrices();
        double[][] complete = informationMatrices[0];
        double[][] jMat = informationMatrices[1];
        // precision of P_G(\beta | g, Z, \tau)
        double[][] precisionMarginal = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, complete);
        // precision of P_G(\beta | \gamma, g, Z, \tau)
        double[][] precisionConditional = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, jMat);

        double logDetMarginal = GMRFDenseMatrixUtils.logDetSPD(precisionMarginal);
        double logDetConditional = GMRFDenseMatrixUtils.logDetSPD(precisionConditional);
        double mutualInf = 0.5*(logDetMarginal - logDetConditional);

        if(mutualInf < -1e-6){
            throw  new RuntimeException("Approximate mutual information of log effective population" +
                    "size and effect size coefficients (conditional on genealogy, covariates and precision)" +
                    "is negative. This should not happen.");
        }
        return Math.max(mutualInf, 0);
    }

    // Computes \mu_{\beta | \gamma, g, Z, \tau}
    private double[] computeBetaGivenGammaMean(){
        double[][] complete = computeInformationMatrices()[0];
        double[][] precision = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, complete);
        double[] tauZQGamma = tauZTransposeQTimes(gammaAsArray());
        double[] betaPriorPrecTimesMean = GMRFDenseMatrixUtils.multMatVec(betaPriorPrecision,betaPriorMean);
        double[] tempVect = new double[betaPriorPrecTimesMean.length];
        for (int j = 0; j < tempVect.length; j++) {
            tempVect[j] = betaPriorPrecTimesMean[j] + tauZQGamma[j];
        }
        return GMRFDenseMatrixUtils.solveSPD(precision, tempVect);
    }

    // Computes diagonal of \Sigma_{\beta | \gamma, g, Z, \tau}
    private double[] computeBetaGivenGammaVarianceDiag(){
        double[][] complete = computeInformationMatrices()[0];
        double[][] precision = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, complete);
        return GMRFDenseMatrixUtils.diagonalOfInverseSPD(precision);
    }

    private double[] computeBetaTildeMean(){
        double[][] jMat = computeInformationMatrices()[1];
        double[][] precision = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, jMat);

        double[] sigmaGInverseMuG = approximation.getSigmaGInverseMuG();
        double[] u = approximation.solveSystemVector(sigmaGInverseMuG);
        double[] tauZQu = tauZTransposeQTimes(u);
        double[] betaPriorPrecTimesMean = GMRFDenseMatrixUtils.multMatVec(betaPriorPrecision,betaPriorMean);

        double[] tempVect = new double[betaPriorPrecTimesMean.length];
        for (int j = 0; j < tempVect.length; j++) {
            tempVect[j] = betaPriorPrecTimesMean[j] + tauZQu[j];
        }
        return GMRFDenseMatrixUtils.solveSPD(precision, tempVect);
    }


    private double[] computeBetaTildeVarianceDiag(){
        double[][] jMat = computeInformationMatrices()[1];
        double[][] precision = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, jMat);
        return GMRFDenseMatrixUtils.diagonalOfInverseSPD(precision);
    }

    private double computeTauConditionalShape(){
        int m = approximation.getGammaHat().length-1;
        return tauShape+m/2.0;
    }

    private double computeTauConditionalRate(){
        double[] gamma = gammaAsArray();
        double[] zBeta = approximation.getLikelihood().getZBetaVect();
        double[] diff = new double[gamma.length];
        for (int i = 0; i < gamma.length; i++) {
            diff[i] = gamma[i]-zBeta[i];
        }

        SymmTridiagMatrix rawQ = approximation.getRawQ();
        DenseVector diffVect = new DenseVector(diff);
        DenseVector qDiffVect = new DenseVector(diff.length);
        rawQ.mult(diffVect, qDiffVect);
        double quadraticForm = diffVect.dot(qDiffVect);

        return tauRate + 0.5*quadraticForm;
    }

    // Per-interval Lagrange remainder bound on second-order Taylor approximation error
    // |remainder_k| <= (1/6)*SS_k*\exp{-min(\gamma_k, \gammaHat_k)}*|\gamma_k - \gammaHat_k|^3
    private double[] computeLagrangeBoundPerInterval(){
        double[] gammaHat = approximation.getGammaHat();
        double[] suffStat = approximation.getSuffStat();

        double[] bound = new double[gammaHat.length];
        for(int k = 0; k < gammaHat.length; k++){
            double gammaK = gamma.getParameterValue(k);
            double diff = gammaK - gammaHat[k];
            double minVal = Math.min(gammaK, gammaHat[k]);
            bound[k] = (1.0/6.0)*suffStat[k]*Math.exp(-minVal)*Math.pow(Math.abs(diff),3);
        }
        return bound;
    }


    public LogColumn[] getColumns(){
        LogColumn[] columns = new LogColumn[getDimension()];
        for (int i = 0; i < columns.length; i++) {
            int dim = i;
            String label = (getDimension() == 1) ? getStatisticName() : getStatisticName() + (dim +1);
            columns[i] = new NumberColumn(label) {
                @Override
                public double getDoubleValue() {
                    return getStatisticValue(dim);
                }
            };
        }
        return columns;
    }

}