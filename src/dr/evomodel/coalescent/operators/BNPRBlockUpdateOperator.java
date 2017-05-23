package dr.evomodel.coalescent.operators;

import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.BNPRLikelihood;
import dr.evomodelxml.coalescent.operators.BNPRBlockUpdateOperatorParser;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import no.uib.cipr.matrix.*;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRBlockUpdateOperator extends SimpleMCMCOperator { // implements GibbsOperator

    private int fieldLength;

    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private List<Parameter> betaParameter;
    private List<MatrixParameter> covariates;

    BNPRLikelihood bnprField;

    private double[] zeros;

    public static final double TWO_TIMES_PI =6.283185;

    public BNPRBlockUpdateOperator(BNPRLikelihood bnprLikelihood, double weight) {
        bnprField = bnprLikelihood;
        popSizeParameter = bnprLikelihood.getPopSizeParameter();
        precisionParameter = bnprLikelihood.getPrecisionParameter();

//        betaParameter = bnprLikelihood.getBetaListParameter();
//        covariates = bnprLikelihood.getCovariates();

        fieldLength = popSizeParameter.getDimension();

        setWeight(weight);

        this.bnprField = bnprLikelihood;

        this.zeros = new double[fieldLength];
    }

    public DenseVector ESS() {
        DenseVector currentPopSize = new DenseVector(bnprField.getPopSizeParameter().getParameterValues());
        double currentLoglik = bnprField.getLogLikelihood();

        // Choose ellipse
        SymmTridiagMatrix currentQ = bnprField.getScaledWeightMatrix(precisionParameter.getParameterValue(0));
        BandCholesky chol = BandCholesky.factorize(new UpperSPDBandMatrix(currentQ, 1));
        DenseVector nu = getMultiNormal(new DenseVector(zeros), chol);
//        DenseVector nu = new DenseVector(fieldLength);

        // Log-likelihood threshold
        double u = MathUtils.nextDouble();
        double logy = currentLoglik + Math.log(u);

        // Draw a initial proposal, also defining a bracket
        double t = TWO_TIMES_PI + MathUtils.nextDouble();
        double tMin = t - TWO_TIMES_PI;
        double tMax = t;

        DenseVector q = currentPopSize.copy();
        q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

        double l = bnprField.getLogLikelihoodSubGamma(q.getData());

        while (l < logy) {

            if (t < 0) {
                tMin = t;
            } else {
                tMax = t;
            }

            t = tMin + MathUtils.nextDouble() * (tMax - tMin);
            q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

            l = bnprField.getLogLikelihoodSubGamma(q.getData());
        }

        return q;
    }

    public double getNewPrecision(double currentValue) {
        //double length = scaleFactor - 1 / scaleFactor;
        double returnValue;
        double shape, scale;
        DenseVector currentGamma = new DenseVector(bnprField.getPopSizeParameter().getParameterValues());

        shape = bnprField.getPrecAlpha() + (bnprField.getNumGridPoints() - 1) / 2.0; // TODO: Watch for off-by-one
        scale = bnprField.getPrecBeta() + currentGamma.dot(bnprField.getScaledWeightMatrix(currentValue).mult(currentGamma, new DenseVector(currentGamma.size())));

        returnValue = GammaDistribution.nextGamma(shape, scale);

        return returnValue;
    }

    public DenseVector getMultiNormalMean(DenseVector canonVector, BandCholesky cholesky) {

        DenseVector tempValue = new DenseVector(zeros);
        DenseVector mean = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = cholesky.getU();

        // Assume Cholesky factorization of the precision matrix Q = LL^T

        // 1. Solve L\omega = b

        CholeskyUpper.transSolve(canonVector, tempValue);

        // 2. Solve L^T \mu = \omega

        CholeskyUpper.solve(tempValue, mean);

        return mean;
    }

    public DenseVector getMultiNormal(DenseVector standNorm, DenseVector mean, BandCholesky cholesky) {

        DenseVector returnValue = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = cholesky.getU();

        // 3. Solve L^T v = z

        CholeskyUpper.solve(standNorm, returnValue);

        // 4. Return x = \mu + v

        returnValue.add(mean);

        return returnValue;
    }

    public static DenseVector getMultiNormal(DenseVector Mean, BandCholesky cholesky) {
        int length = Mean.size();
        DenseVector tempValue = new DenseVector(length);

        for (int i = 0; i < length; i++)
            tempValue.set(i, MathUtils.nextGaussian());

        DenseVector returnValue = new DenseVector(Mean.size());
        UpperTriangBandMatrix CholeskyUpper = cholesky.getU();

        // 3. Solve L^T v = z
        CholeskyUpper.solve(tempValue, returnValue);
        // 4. Return x = \mu + v
        returnValue.add(Mean);
        return returnValue;
    }

    public static DenseVector getMultiNormal(DenseVector mean, UpperSPDDenseMatrix variance) {
        int length = mean.size();
        DenseVector tempValue = new DenseVector(length);
        DenseVector returnValue = new DenseVector(length);
        UpperSPDDenseMatrix ab = variance.copy();

        for (int i = 0; i < returnValue.size(); i++)
            tempValue.set(i, MathUtils.nextGaussian());

        DenseCholesky chol = new DenseCholesky(length, true);
        chol.factor(ab);

        UpperTriangDenseMatrix x = chol.getU();

        x.transMult(tempValue, returnValue);
        returnValue.add(mean);
        return returnValue;
    }

    public double doOperation() {
        double currentPrecision = precisionParameter.getParameterValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision);

        //double currentLambda = this.lambdaParameter.getParameterValue(0);
        //double proposedLambda = this.getNewLambda(currentLambda, lambdaScaleFactor);

        precisionParameter.setParameterValue(0, proposedPrecision);
        //lambdaParameter.setParameterValue(0, proposedLambda);

        DenseVector currentGamma = new DenseVector(bnprField.getPopSizeParameter().getParameterValues());
        DenseVector proposedGamma = new DenseVector(ESS());

        for (int i = 0; i < fieldLength; i++)
            popSizeParameter.setParameterValueQuietly(i, proposedGamma.get(i));

        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();

        double hRatio = 0; // This step always accepts.

        return hRatio;
    }

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return BNPRBlockUpdateOperatorParser.BNPR_BLOCK_OPERATOR;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {
        return "This operator should not need tuning, and should accept with probability 1.";
    }

    // Main function for testing

    public static void main(String[] args) throws Exception {
//        NewickImporter importer = new NewickImporter("((((5:0.5,1:0.2):0.5,0:1):0.2,2:0.8):0.2,3:1.4)");
//        System.out.println(System.getProperty("user.home"));
        NexusImporter importer = new NexusImporter(new FileReader(System.getProperty("user.home") + "/Dropbox/VladResearch/Julia/influenza_max_clade.nex"));
        Tree tree = importer.importNextTree();

        System.out.println("Got here A");

        System.out.println("Number of nodes: " + tree.getNodeCount());

        BNPRLikelihood bnprLikelihood = new BNPRLikelihood(tree, new Parameter.Default(1.0),
                null, new Parameter.Default(1.0), 679.66184, 100, false);

        System.out.println("Got here B");

        BNPRBlockUpdateOperator bnprBlockUpdateOperator = new BNPRBlockUpdateOperator(bnprLikelihood, 1.0);

        System.out.println("Got here C");

        System.out.println("Pop sizes: " + Arrays.toString(bnprBlockUpdateOperator.bnprField.getPopSizeParameter().getParameterValues()));

        System.out.println("Got here D");

        bnprBlockUpdateOperator.doOperation();
        System.out.println("Pop sizes: " + Arrays.toString(bnprBlockUpdateOperator.bnprField.getPopSizeParameter().getParameterValues()));

        System.out.println("Got here E");

        for (int i = 0; i < 10; i++) {
            bnprBlockUpdateOperator.doOperation();
        }

        System.out.println("Opening file for writing...");
        PrintWriter pw = new PrintWriter(new FileWriter(System.getProperty("user.home") + "/Documents/Research/BEASText/doOperationOutput.txt"));

        for (int i = 0; i < 30; i++) {
            bnprBlockUpdateOperator.doOperation();
            if (i % 10 == 9) {
                System.out.println("Iteration " + i);
                pw.println(Arrays.toString(bnprBlockUpdateOperator.bnprField.getPopSizeParameter().getParameterValues()));
            }
        }

        pw.close();
    }
}
