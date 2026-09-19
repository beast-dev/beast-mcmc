package dr.evomodel.coalescent;

import dr.xml.*;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.util.FileHelpers;
import org.apache.commons.math.special.Gamma;

import java.io.*;
import java.nio.file.Paths;
import java.util.List;

public class SkygridMutualInformationAnalysis {

    public enum Mode {
        FIXED_TAU,
        CONDITIONAL_TAU_IS,
        MARGINAL_TAU_AVERAGE,
        FULLY_MARGINAL
    }

    private final Mode mode;
    private final double[][] gammaSamples;
    private final double[][] betaSamples;
    private final double[] tauSamples;
    private final double fixedTauValue;

    private final double[][] zMat;
    private final double[][] betaPriorPrecision;
    private final double[] betaPriorMean;
    private final double tauShape;
    private final double tauRate;
    private final double[] tau0List;

    private final int numQuadraturePoints;
    private final double logTauLower;
    private final double logTauUpper;

    private final int N;
    private final int M; // numGridPoints (fieldLength - 1)
    private final int P;

    private final double[] gammaQgamma;
    private final double[][] zqGamma;
    private final double[][] zqz;
    private final double[] bMat0invb0Vec;
    private final double b0bMat0invb0;

    private String report;

    public SkygridMutualInformationAnalysis(Mode mode,
                                            double[][] gammaSamples,
                                            double[][] betaSamples,
                                            double[] tauSamples,
                                            double fixedTauValue,
                                            double[][] zMat,
                                            double[][] rawQ,
                                            double[][] betaPriorPrecision,
                                            double[] betaPriorMean,
                                            double tauShape,
                                            double tauRate,
                                            double[] tau0List,
                                            int numQuadraturePoints,
                                            double logTauLower,
                                            double logTauUpper){
        this.mode = mode;
        this.gammaSamples = gammaSamples;
        this.betaSamples = betaSamples;
        this.tauSamples = tauSamples;
        this.fixedTauValue = fixedTauValue;
        this.zMat = zMat;
        this.betaPriorPrecision = betaPriorPrecision;
        this.betaPriorMean = betaPriorMean;
        this.tauShape = tauShape;
        this.tauRate = tauRate;
        this.tau0List = tau0List;
        this.numQuadraturePoints = numQuadraturePoints;
        this.logTauLower = logTauLower;
        this.logTauUpper = logTauUpper;

        this.N = gammaSamples.length;
        this.M = zMat.length-1;
        this.P = zMat[0].length;

        // compute Z'QZ and, for eeach sample, compute \gamma'Q\gamma Z'Q\gamma

        this.zqz = computeZTransposeQZ(rawQ, zMat);
        this.gammaQgamma = new double[N];
        this.zqGamma = new double[N][P];
        for(int i = 0; i < N; i++){
            double[] qGamma = multiplyByQ(rawQ, gammaSamples[i]);
            double gqg = 0.0;
            for (int k = 0; k < (M+1); k++){
                gqg += gammaSamples[i][k]*qGamma[k];
            }
            gammaQgamma[i] = gqg;
            zqGamma[i] = GMRFDenseMatrixUtils.transposeMultiplyVector(zMat, qGamma);
        }

        this.bMat0invb0Vec = GMRFDenseMatrixUtils.multMatVec(betaPriorPrecision, betaPriorMean);
        double sum = 0.0;
        for(int j = 0; j < P; j++){
            sum += betaPriorMean[j]*bMat0invb0Vec[j];
        }
        this.b0bMat0invb0 = sum;
    }

    private static double[] multiplyByQ(double[][] rawQ, double[] v){
        int n = v.length;
        double[] result = new double[n];
        for(int i = 0; i < n; i++){
            double s = 0.0;
            for(int j = 0; j < n; j++){
                s += rawQ[i][j]*v[j];
            }
            result[i] = s;
        }
        return result;
    }

    private static double[][] computeZTransposeQZ(double[][] rawQ, double[][] zMat){
        int n = zMat.length;
        int p = zMat[0].length;
        double[][] qz = new double[n][p];
        for(int j = 0; j < p; j++){
            double[] col = new double[n];
            for(int i = 0; i < n; i++){
                col[i] = zMat[i][j];
            }
            double[] qCol = multiplyByQ(rawQ, col);
            for(int i = 0; i < n; i++){
                qz[i][j] = qCol[i];
            }
        }
        return GMRFDenseMatrixUtils.transposeMultiply(zMat, qz);
    }

    private double quadraticFormGammaBeta(int i){
        double[] beta = betaSamples[i];
        double cross = 0.0;
        for (int j = 0; j < P; j++){
            cross += beta[j]*zqGamma[i][j];
        }
        double betaZQZBeta = GMRFDenseMatrixUtils.quadraticForm(zqz, beta);
        return gammaQgamma[i] - 2.0*cross + betaZQZBeta;
    }

    private double logGammaDensity(double x, double shape, double rate){
        return shape*Math.log(rate) - Gamma.logGamma(shape) +
                (shape-1.0)*Math.log(x)-rate*x;
    }

    private static double logSumExp(double[] logValues){
        double max = Double.NEGATIVE_INFINITY;
        for(double v : logValues){
            if(v > max){
                max = v;
            }
        }
        if(Double.isInfinite(max)){
            return max;
        }
        double sum = 0.0;
        for (double v : logValues){
            sum += Math.exp(v-max);
        }
        return max + Math.log(sum);
    }

    // normalize the log weights into weights that sum to 1
    private static double[] normalizeLogWeights(double[] logWeights){
        double lse = logSumExp(logWeights);
        double[] w = new double[logWeights.length];
        for(int i = 0; i < w.length; i++){
            w[i] = Math.exp(logWeights[i]-lse);
        }
        return w;
    }

    private static double effectiveSampleSize(double[] normalizedWeights){
        double sumSq = 0.0;
        for (double w : normalizedWeights){
            sumSq += w*w;
        }
        return 1.0/sumSq;
    }

    private double[] rhsForTau(int sampleIndex, double tau){
        double[] rhs = new double[P];
        for (int j = 0; j < P; j++){
            rhs[j] = bMat0invb0Vec[j] + tau*zqGamma[sampleIndex][j];
        }
        return rhs;
    }

    private double closedFormConditionalEntropy(double logDetPrecision){
        return 0.5*P*Math.log(2*Math.PI*Math.E)-0.5*logDetPrecision;
    }

    // For a specific \tau value, has the precision B_0^{-1}\tau*Z'QZ,
    // its Cholesky factor, and log determinant
    private static class BetaGivenGammaSystem{
        final double[][] prec;
        final double[][] cholesky;
        final double logDet;

        BetaGivenGammaSystem(double[][] betaPriorPrecision, double[][] zqz, double tau){
            prec = GMRFDenseMatrixUtils.addMat(betaPriorPrecision, GMRFDenseMatrixUtils.scaleMult(zqz, tau));
            cholesky = GMRFDenseMatrixUtils.choleskyLower(prec);
            double logDetAcc = 0.0;
            for (int i = 0; i < cholesky.length; i++){
                logDetAcc += Math.log(cholesky[i][i]);
            }
            logDet = 2.0*logDetAcc;
        }

        double[] mean(double[] rhs){
            return GMRFDenseMatrixUtils.solveGivenCholeskyLower(cholesky, rhs);
        }
    }


    // Estimating the mutual information for a fixed \tau
    private double[] runFixedTau(){
        double tau = fixedTauValue;
        if(tauSamples != null){
            double minTau = tauSamples[0];
            double maxTau = tauSamples[0];
            for (double t : tauSamples){
                minTau = Math.min(minTau, t);
                maxTau = Math.max(maxTau, t);
            }
            if(maxTau - minTau > 1e-6*(1.0 + Math.abs(maxTau))){
                throw new RuntimeException("We are under fixed tau analyis, but tau is not being held constant,");
            }
            tau = tauSamples[0];
        }

        BetaGivenGammaSystem sys = new BetaGivenGammaSystem(betaPriorPrecision, zqz, tau);
        double hCond = closedFormConditionalEntropy(sys.logDet);

        double[][] means = new double[N][];
        for (int s = 0; s < N; s++){
            means[s] = sys.mean(rhsForTau(s, tau));
        }

        double sumLogMixture = 0.0;
        for(int r = 0; r < N; r++){
            double[] logDensities = new double[N];
            for(int s = 0; s < N; s++){
                logDensities[s] = GMRFDenseMatrixUtils.logMvnDensityGivenPrecision(
                        betaSamples[r], means[s], sys.prec, sys.logDet);

            }
            sumLogMixture += logSumExp(logDensities)-Math.log(N);
        }
        double hMarg = -sumLogMixture/N;

        return new double[]{hMarg-hCond, tau, Double.NaN};
    }

    // conditional \tau importance sampling
    private double[] estimateAtTau0(double tau0){
        double[] logWeights = new double[N];
        for (int r =0; r < N; r++){
            double rate = tauRate + 0.5*quadraticFormGammaBeta(r);
            logWeights[r] = logGammaDensity(tau0, tauShape + M/2.0, rate);
        }
        double[] w = normalizeLogWeights(logWeights);
        double ess = effectiveSampleSize(w);

        BetaGivenGammaSystem sys = new BetaGivenGammaSystem(betaPriorPrecision, zqz, tau0);
        double hCond = closedFormConditionalEntropy(sys.logDet);

        double[][] means = new double[N][];
        for(int r = 0; r < N; r++){
            means[r] = sys.mean(rhsForTau(r, tau0));
        }

        double hMarg = 0.0;
        for(int r = 0; r < N; r++){
            double[] logDensities = new double[N];
            for(int rPrime = 0; rPrime < N; rPrime++){
                logDensities[rPrime] = Math.log(w[rPrime])
                        + GMRFDenseMatrixUtils.logMvnDensityGivenPrecision(
                                betaSamples[r], means[rPrime], sys.prec, sys.logDet);
            }
            double logMixtureR = logSumExp(logDensities);
            hMarg += -w[r]*logMixtureR;
        }
        return new double[]{hMarg-hCond, tau0, ess};
    }


    // marginal \tau average
    private double[] runMarginalTauAverage(){
        double[] rates = new double[N];
        for (int r = 0; r < N; r++){
            rates[r] = tauRate + 0.5*quadraticFormGammaBeta(r);
        }

        double sumHCond = 0.0;
        double sumHMargContribution = 0.0;
        double sumEss = 0.0;

        for(int s = 0; s < N; s++){
            double tauS = tauSamples[s];
            BetaGivenGammaSystem sys = new BetaGivenGammaSystem(betaPriorPrecision, zqz, tauS);
            sumHCond += closedFormConditionalEntropy(sys.logDet);

            double[] logWeights = new double[N];
            for (int r = 0; r < N; r++){
                logWeights[r] = logGammaDensity(tauS, tauShape + M/2.0, rates[r]);
            }
            double[] w = normalizeLogWeights(logWeights);
            sumEss += effectiveSampleSize(w);

            double[] logDensities = new double[N];
            for (int r = 0; r < N; r++){
                double[] meanR = sys.mean(rhsForTau(r, tauS));
                logDensities[r] = Math.log(w[r]) + GMRFDenseMatrixUtils.logMvnDensityGivenPrecision(
                    betaSamples[s], meanR, sys.prec, sys.logDet);
            }
            double logPBetaSGivenTauS = logSumExp(logDensities);
            sumHMargContribution += -logPBetaSGivenTauS;
        }

        double hCondAvg = sumHCond/N;
        double hMargAvg = sumHMargContribution/N;
        double essAvg = sumEss/N;

        return new double[]{hMargAvg-hCondAvg, Double.NaN, essAvg};
    }

    // fully marginal mutual information

    private static double[] subtract(double[] a, double[] b){
        double[] result = new double[a.length];
        for(int i =0; i < a.length; i++){
            result[i] = a[i] - b[i];
        }
        return result;
    }

    // matching nodes and weights for Simpson's rule
    // \int f \approx \sum_k weights[k]*f(points[k])
    // numPoints is shifted up by one if it is even, since we
    // need an even number of intervals (and odd number points)
    private static double[][] simpsonNodesAndWeights(int numPoints, double lower, double upper) {
        if (numPoints % 2 == 0) {
            numPoints += 1;
        }
        int n = numPoints - 1; // number of intervals, even
        double h = (upper-lower)/n;

        double[] points = new double[numPoints];
        double[] weights = new double[numPoints];
        for (int i = 0; i <= n; i++) {
            points[i] = lower + i*h;
        }
        weights[0] = h/3.0;
        weights[n] = h/3.0;
        for (int i = 1; i < n; i++) {
            weights[i] = (i%2 == 1) ? (4.0*h/3.0):(2.0*h/3.0);
        }
        return new double[][]{points, weights};
    }

    private static class QuadratureNode {
        final double u, tau, logWeight, logDetM;
        final double[][] choleskyM;

        QuadratureNode(double u, double weight, double[][] betaPriorPrecision, double[][] zqz) {
            this.u = u;
            this.tau = Math.exp(u);
            this.logWeight = Math.log(weight);
            double[][] mMat = GMRFDenseMatrixUtils.addMat(betaPriorPrecision,
                    GMRFDenseMatrixUtils.scaleMult(zqz, tau));
            this.choleskyM = GMRFDenseMatrixUtils.choleskyLower(mMat);
            double acc = 0.0;
            for (int i = 0; i < choleskyM.length; i++) {
                acc += Math.log(choleskyM[i][i]);
            }
            this.logDetM = 2.0*acc;
        }
    }

    private double runFullyMarginal() {
        double[][] nodesAndWeights = simpsonNodesAndWeights(numQuadraturePoints, logTauLower, logTauUpper);
        double[] points = nodesAndWeights[0];
        double[] weights = nodesAndWeights[1];

        QuadratureNode[] nodes = new QuadratureNode[points.length];
        for (int k = 0; k < points.length; k++) {
            nodes[k] = new QuadratureNode(points[k], weights[k], betaPriorPrecision, zqz);
        }

        double[] logPGammaGivenZ = new double[N];
        for (int i = 0; i < N; i++) {
            double[] logIntegrand = new double[nodes.length];
            for (int k = 0; k < nodes.length; k++) {
                QuadratureNode node = nodes[k];
                double[] rhs = rhsForTau(i, node.tau);
                double[] solved = GMRFDenseMatrixUtils.solveGivenCholeskyLower(node.choleskyM, rhs);
                double quad = 0.0;
                for (int j = 0; j < P; j++) {
                    quad += rhs[j]*solved[j];
                }
                double exponent = -0.5*b0bMat0invb0 - 0.5*node.tau*gammaQgamma[i] + 0.5*quad;
                logIntegrand[k] = (tauShape + M / 2.0)*node.u - tauRate*node.tau
                        - 0.5*node.logDetM + exponent + node.logWeight;
            }
            logPGammaGivenZ[i] = logSumExp(logIntegrand);
        }

        BetaGivenGammaSystem[] systems = new BetaGivenGammaSystem[N];
        double[][] means = new double[N][];
        for (int s = 0; s < N; s++) {
            systems[s] = new BetaGivenGammaSystem(betaPriorPrecision, zqz, tauSamples[s]);
            means[s] = systems[s].mean(rhsForTau(s, tauSamples[s]));
        }

        double logGammaHalfMPlusA = Gamma.logGamma(tauShape + M/2.0);

        double sumLogRatio = 0.0;
        for (int i = 0; i < N; i++) {
            double qi = quadraticFormGammaBeta(i);
            double logNumerator = -0.5*P*Math.log(2*Math.PI)
                    - 0.5 * GMRFDenseMatrixUtils.quadraticForm(betaPriorPrecision, subtract(betaSamples[i], betaPriorMean))
                    + logGammaHalfMPlusA
                    - (tauShape + M / 2.0)*Math.log(tauRate + 0.5*qi)
                    - logPGammaGivenZ[i];

            double[] logDensities = new double[N];
            for (int s = 0; s < N; s++) {
                logDensities[s] = GMRFDenseMatrixUtils.logMvnDensityGivenPrecision(
                        betaSamples[i], means[s], systems[s].prec, systems[s].logDet);
            }
            double logDenominator = logSumExp(logDensities)-Math.log(N);

            sumLogRatio += (logNumerator-logDenominator);
        }
        return sumLogRatio/N;
    }


    // report

    public String getReport(){
        if (report == null){
            run();
        }
        return report;
    }

    private void run() {
        StringBuilder sb = new StringBuilder();
        sb.append("Skygrid-GLM mutual information analysis (mode=").append(mode).append(")\n");
        sb.append("N = ").append(N).append(" retained samples, P = ").append(P)
                .append(" covariates, M = ").append(M).append(" grid points\n\n");

        switch (mode) {
            case FIXED_TAU: {
                double[] result = runFixedTau();
                sb.append("tau (fixed) = ").append(result[1]).append("\n");
                sb.append("I(gamma,beta|g,Z,tau) = ").append(result[0]).append("\n");
                break;
            }
            case CONDITIONAL_TAU_IS: {
                sb.append("tau0\tI(gamma,beta|g,Z,tau0)\tESS\n");
                for (double tau0 : tau0List) {
                    double[] result = estimateAtTau0(tau0);
                    sb.append(result[1]).append("\t").append(result[0]).append("\t").append(result[2]).append("\n");
                }
                break;
            }
            case MARGINAL_TAU_AVERAGE: {
                double[] result = runMarginalTauAverage();
                sb.append("E_tau[I(gamma,beta|g,Z,tau)] = ").append(result[0]).append("\n");
                sb.append("mean importance-sampling ESS across rows = ").append(result[2]).append("\n");
                break;
            }
            case FULLY_MARGINAL: {
                double value = runFullyMarginal();
                sb.append("I(gamma,beta|g,Z) = ").append(value).append("\n");
                break;
            }
        }

        report = sb.toString();
    }

    // parser
    public static final String SKYGRID_MI_ANALYSIS = "skygridMutualInformationAnalysis";
    public static final String MODE = "mode";
    public static final String RESULT_FILE_NAME = "resultsFileName";
    public static final String BURN_IN = "burnIn";
    public static final String BETA_PRIOR_PRECISION = "betaPriorPrecision";
    public static final String BETA_PRIOR_MEAN = "betaPriorMean";
    public static final String TAU_SHAPE = "tauShape";
    public static final String TAU_RATE = "tauRate";
    public static final String FIXED_TAU_VALUE = "fixedTauValue";
    public static final String TAU0_LIST = "tau0List";
    public static final String NUM_QUADRATURE_POINTS = "numQuadraturePoints";
    public static final String LOG_TAU_LOWER = "logTauLower";
    public static final String LOG_TAU_UPPER = "logTauUpper";

    private static final int DEFAULT_NUM_QUADRATURE_POINTS = 401;
    private static final double DEFAULT_LOG_TAU_LOWER = -15.0;
    private static final double DEFAULT_LOG_TAU_UPPER = 15.0;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SKYGRID_MI_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            GMRFMultilocusSkyrideLikelihood likelihood =
                    (GMRFMultilocusSkyrideLikelihood) xo.getChild(GMRFMultilocusSkyrideLikelihood.class);
            if (likelihood == null) {
                throw new XMLParseException(SKYGRID_MI_ANALYSIS + " requires a GMRFMultilocusSkyrideLikelihood.");
            }

            String modeString = xo.getStringAttribute(MODE);
            SkygridMutualInformationAnalysis.Mode mode;
            switch (modeString) {
                case "fixedTau":
                    mode = Mode.FIXED_TAU;
                    break;
                case "conditionalTauIS":
                    mode = Mode.CONDITIONAL_TAU_IS;
                    break;
                case "marginalTauAverage":
                    mode = Mode.MARGINAL_TAU_AVERAGE;
                    break;
                case "fullyMarginal":
                    mode = Mode.FULLY_MARGINAL;
                    break;
                default:
                    throw new XMLParseException("Unrecognized mode: \"" + modeString + "\". Expected one of: " +
                            "fixedTau, conditionalTauIS, marginalTauAverage, fullyMarginal.");
            }

            double[][] betaPriorPrecision;
            {
                XMLObject cxo = xo.getChild(BETA_PRIOR_PRECISION);
                dr.inference.model.MatrixParameter mp =
                        (dr.inference.model.MatrixParameter) cxo.getChild(dr.inference.model.MatrixParameter.class);
                int rows = mp.getRowDimension();
                int cols = mp.getColumnDimension();
                betaPriorPrecision = new double[rows][cols];
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        betaPriorPrecision[i][j] = mp.getParameterValue(i, j);
                    }
                }
            }

            double[] betaPriorMean;
            {
                XMLObject cxo = xo.getChild(BETA_PRIOR_MEAN);
                dr.inference.model.Parameter p =
                        (dr.inference.model.Parameter) cxo.getChild(dr.inference.model.Parameter.class);
                betaPriorMean = p.getParameterValues();
            }

            double tauShape;
            {
                XMLObject cxo = xo.getChild(TAU_SHAPE);
                dr.inference.model.Parameter p =
                        (dr.inference.model.Parameter) cxo.getChild(dr.inference.model.Parameter.class);
                tauShape = p.getParameterValue(0);
            }

            double tauRate;
            {
                XMLObject cxo = xo.getChild(TAU_RATE);
                dr.inference.model.Parameter p =
                        (dr.inference.model.Parameter) cxo.getChild(dr.inference.model.Parameter.class);
                tauRate = p.getParameterValue(0);
            }

            double fixedTauValue = xo.getAttribute(FIXED_TAU_VALUE, Double.NaN);
            if (mode == Mode.FIXED_TAU && Double.isNaN(fixedTauValue)) {
                throw new XMLParseException("mode=\"fixedTau\" requires a fixedTauValue attribute, used if no " +
                        "tau column is present in the log file (or as a cross-check if one is).");
            }

            double[] tau0List = null;
            if (mode == Mode.CONDITIONAL_TAU_IS) {
                if (!xo.hasChildNamed(TAU0_LIST)) {
                    throw new XMLParseException("mode=\"conditionalTauIS\" requires a <tau0List> child element.");
                }
                XMLObject cxo = xo.getChild(TAU0_LIST);
                dr.inference.model.Parameter p =
                        (dr.inference.model.Parameter) cxo.getChild(dr.inference.model.Parameter.class);
                tau0List = p.getParameterValues();
            }

            int numQuadraturePoints = xo.getAttribute(NUM_QUADRATURE_POINTS, DEFAULT_NUM_QUADRATURE_POINTS);
            double logTauLower = xo.getAttribute(LOG_TAU_LOWER, DEFAULT_LOG_TAU_LOWER);
            double logTauUpper = xo.getAttribute(LOG_TAU_UPPER, DEFAULT_LOG_TAU_UPPER);

            String fileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
            String resultFileName = xo.hasAttribute(RESULT_FILE_NAME) ? xo.getStringAttribute(RESULT_FILE_NAME) : null;
            int burnIn = xo.getAttribute(BURN_IN, 0);

            try {
                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();
                if (!file.isAbsolute()) {
                    parent = (parent == null) ? System.getProperty("user.dir")
                            : Paths.get(System.getProperty("user.dir"), parent).toString();
                }
                file = new File(parent, name);

                LogFileTraces traces = new LogFileTraces(file.getAbsolutePath(), file);
                traces.loadTraces();
                traces.setBurnIn(burnIn);

                dr.inference.model.Parameter popSizeParameter = likelihood.getPopSizeParameter();
                int fieldLength = popSizeParameter.getDimension();
                double[][] gammaSamples = readMultiDimensional(traces, popSizeParameter.getId(), fieldLength);

                List<dr.inference.model.Parameter> betaList = likelihood.getBetaListParameter();
                int pDim = 0;
                for (dr.inference.model.Parameter b : betaList) {
                    pDim += b.getDimension();
                }
                double[][] betaSamples = readBetaList(traces, betaList);

                double[] tauSamples = null;
                dr.inference.model.Parameter precisionParameter = likelihood.getPrecisionParameter();
                if (traceExists(traces, precisionParameter.getId())) {
                    tauSamples = readSingleColumn(traces, precisionParameter.getId());
                } else if (mode != Mode.FIXED_TAU) {
                    throw new XMLParseException("mode=\"" + modeString + "\" requires a logged tau column " +
                            "(\"" + precisionParameter.getId() + "\"), which was not found in the log file.");
                }

                double[][] zMat = likelihood.getDesignMat();
                double lambda = likelihood.getLambdaParameter().getParameterValue(0);
                double[][] rawQ = toDoubleArray(likelihood.getScaledWeightMatrix(1.0, lambda), fieldLength);

                SkygridMutualInformationAnalysis analysis = new SkygridMutualInformationAnalysis(
                        mode, gammaSamples, betaSamples, tauSamples, fixedTauValue,
                        zMat, rawQ, betaPriorPrecision, betaPriorMean, tauShape, tauRate,
                        tau0List, numQuadraturePoints, logTauLower, logTauUpper);

                String reportText = analysis.getReport();
                System.out.println(reportText);

                if (resultFileName != null) {
                    FileWriter fw = new FileWriter(resultFileName, true);
                    BufferedWriter bw = new BufferedWriter(fw);
                    bw.write(reportText);
                    bw.flush();
                    bw.close();
                }

                return analysis;

            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            } catch (IOException ioe) {
                throw new XMLParseException(ioe.getMessage());
            } catch (TraceException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        private double[][] toDoubleArray(no.uib.cipr.matrix.SymmTridiagMatrix mat, int n) {
            double[][] result = new double[n][n];
            for (int i = 0; i < n; i++) {
                result[i][i] = mat.get(i, i);
                if (i > 0) {
                    result[i][i - 1] = mat.get(i, i - 1);
                    result[i - 1][i] = mat.get(i - 1, i);
                }
            }
            return result;
        }

        private boolean traceExists(LogFileTraces traces, String name) {
            for (int i = 0; i < traces.getTraceCount(); i++) {
                if (traces.getTraceName(i).trim().equals(name)) {
                    return true;
                }
            }
            return false;
        }

        private int findTraceIndex(LogFileTraces traces, String name) throws XMLParseException {
            for (int i = 0; i < traces.getTraceCount(); i++) {
                if (traces.getTraceName(i).trim().equals(name)) {
                    return i;
                }
            }
            throw new XMLParseException("Column '" + name + "' can not be found for " + SKYGRID_MI_ANALYSIS + " element.");
        }

        private double[] readSingleColumn(LogFileTraces traces, String name) throws XMLParseException {
            int idx = findTraceIndex(traces, name);
            List<Double> values = traces.getValues(idx);
            double[] result = new double[values.size()];
            for (int i = 0; i < result.length; i++) {
                result[i] = values.get(i);
            }
            return result;
        }

        private double[][] readMultiDimensional(LogFileTraces traces, String baseName, int dimension)
                throws XMLParseException {
            double[][] columns = new double[dimension][];
            for (int d = 0; d < dimension; d++) {
                String colName = (dimension == 1) ? baseName : baseName + (d + 1);
                columns[d] = readSingleColumn(traces, colName);
            }
            int n = columns[0].length;
            double[][] result = new double[n][dimension];
            for (int i = 0; i < n; i++) {
                for (int d = 0; d < dimension; d++) {
                    result[i][d] = columns[d][i];
                }
            }
            return result;
        }

        private double[][] readBetaList(LogFileTraces traces, List<dr.inference.model.Parameter> betaList)
                throws XMLParseException {
            int totalP = 0;
            for (dr.inference.model.Parameter b : betaList) {
                totalP += b.getDimension();
            }
            double[][] perEntry = new double[totalP][];
            int col = 0;
            for (dr.inference.model.Parameter b : betaList) {
                int dim = b.getDimension();
                for (int d = 0; d < dim; d++) {
                    String colName = (dim == 1) ? b.getId() : b.getId() + (d + 1);
                    perEntry[col] = readSingleColumn(traces, colName);
                    col++;
                }
            }
            int n = perEntry[0].length;
            double[][] result = new double[n][totalP];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < totalP; j++) {
                    result[i][j] = perEntry[j][i];
                }
            }
            return result;
        }

        public String getParserDescription() {
            return "Post-processing analysis (run after completion of MCMC somulation) computing exact, non-Gaussian-" +
                    "approximated mutual information estimators for the skygrid-GLM effective population " +
                    "size trajectory and effect size coefficients.";
        }

        public Class getReturnType() {
            return SkygridMutualInformationAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME, "The BEAST log file to read."),
                new StringAttributeRule(MODE, "fixedTau | conditionalTauIS | marginalTauAverage | fullyMarginal"),
                new StringAttributeRule(RESULT_FILE_NAME, "Output file for the result.", true),
                AttributeRule.newIntegerRule(BURN_IN, true),
                AttributeRule.newDoubleRule(FIXED_TAU_VALUE, true),
                AttributeRule.newIntegerRule(NUM_QUADRATURE_POINTS, true),
                AttributeRule.newDoubleRule(LOG_TAU_LOWER, true),
                AttributeRule.newDoubleRule(LOG_TAU_UPPER, true),
                new ElementRule(GMRFMultilocusSkyrideLikelihood.class),
                new ElementRule(BETA_PRIOR_PRECISION, new XMLSyntaxRule[]{
                        new ElementRule(dr.inference.model.MatrixParameter.class)}),
                new ElementRule(BETA_PRIOR_MEAN, new XMLSyntaxRule[]{
                        new ElementRule(dr.inference.model.Parameter.class)}),
                new ElementRule(TAU_SHAPE, new XMLSyntaxRule[]{
                        new ElementRule(dr.inference.model.Parameter.class)}),
                new ElementRule(TAU_RATE, new XMLSyntaxRule[]{
                        new ElementRule(dr.inference.model.Parameter.class)}),
                new ElementRule(TAU0_LIST, new XMLSyntaxRule[]{
                        new ElementRule(dr.inference.model.Parameter.class)}, true),
        };
    };

}