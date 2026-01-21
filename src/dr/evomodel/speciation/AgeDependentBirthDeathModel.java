package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;

public class AgeDependentBirthDeathModel extends AbstractModelLikelihood {

    private final Tree tree;
    private final String name;
    private final Parameter birthRate;
    private final Parameter deathRate;

    private final int timeSteps;
    private final double h;

    private final double[][] birthRateMat;
    private final double[][] deathRateMat;
    private final double[][] expmR;
    private final double[] p0;
    private final double[] logS;
    private final double[] logBranchDens;

    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown;

    private double logLikelihood;
    private double storedLogLikelihood;

    public AgeDependentBirthDeathModel(String name,
                                       Tree tree,
                                       Parameter birthRate,
                                       Parameter deathRate,
                                       double originTime,
                                       int timeSteps) {
        super(name);

        this.name = name;
        this.tree = tree;

        // Might want to use jagged arrays for the rate matrices
        this.birthRate = birthRate;
        this.birthRateMat = new double[timeSteps+1][timeSteps+1];

        this.deathRate = deathRate;
        this.deathRateMat = new double[timeSteps+1][timeSteps+1];

        this.expmR = new double[timeSteps+1][timeSteps+1];
        this.p0 = new double[timeSteps+1];
        this.logS = new double[timeSteps+1];
        this.logBranchDens = new double[timeSteps+1];

        this.timeSteps = timeSteps;
        
        this.h = originTime / timeSteps;

        addVariable(birthRate);
        addVariable(deathRate);

        fillBirthRates();
        fillDeathRates();
        fillExpmR();
        fillp0();
        fillLogS();

        likelihoodKnown = false;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected void acceptState() {

    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == birthRate) {
            fillBirthRates();
            likelihoodKnown = false;
        } else if (variable == deathRate) {
            fillDeathRates();
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Unknown parameter");
        }
        // Is this what we want to do?
        fillExpmR();
        fillp0();
        fillLogS();
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private void fillBirthRates() {
        for (int i = 0; i <= timeSteps; i++) {
            for (int j = 0; j <= timeSteps - i; j++) {
                birthRateMat[i][j] = birthRate.getParameterValue(0); // Right now constant rates --> Make more interesting
            }
        }
    }

    private void fillDeathRates() {
        for (int i = 0; i <= timeSteps; i++) {
            for (int j = 0; j <= timeSteps - i; j++) {
                deathRateMat[i][j] = deathRate.getParameterValue(0); // Right now constant rates --> Make more interesting
            }
        }
    }
    
    private void fillExpmR() {
        for (double[] row : expmR) {
            java.util.Arrays.fill(row, 0.0);
        }

        // Compute upper triangular cumulative rates (Row-Major)
        for (int i = timeSteps - 1; i >= 0; i--) {
            for (int j = 1; j <= timeSteps - i; j++) {
                double r1 = birthRateMat[i][j] + deathRateMat[i][j];
                double r2 = birthRateMat[i + 1][j - 1] + deathRateMat[i + 1][j - 1];
                expmR[i][j] = h * (r1 + r2) / 2.0 + expmR[i + 1][j - 1];
            }
        }

        // Exponentiate rates (Row-Major)
        for (int i = 0; i <= timeSteps; i++) {
            for (int j = 0; j <= timeSteps - i; j++) {
                expmR[i][j] = Math.exp(-expmR[i][j]);
            }
        }
    }

    private void fillp0() { // Give a better name maybe? fillExtProb?
        java.util.Arrays.fill(p0, 0.0);

        // Solve integral equation with trapezoidal rule - corresponds to solving a system of quadratic equations
        for (int m = 1; m <= timeSteps; m++) {
            double trap_sum_mu = 0.0;
            double trap_sum_lam = 0.0;

            trap_sum_mu += deathRateMat[0][m] * expmR[0][m] / 2.0;

            for (int i = 1; i < m; i++) {
                int j = m - i;
                trap_sum_mu += deathRateMat[i][j] * expmR[i][j];
                trap_sum_lam += birthRateMat[i][j] * expmR[i][j] * p0[i] * p0[i];
            }

            trap_sum_mu += deathRateMat[m][0] * expmR[m][0] / 2.0;

            trap_sum_mu *= h;
            trap_sum_lam *= h;

            double a = h * birthRateMat[m][0] * expmR[m][0] / 2.0;
            double c = trap_sum_mu + trap_sum_lam;

            if (a == 0.0) {
                p0[m] = c;
            } else {
                double det = 1.0 - 4.0 * a * c;
                if (det < 0.0) {
                    p0[m] = 1.0 / (2.0 * a);
                } else {
                    p0[m] = (1.0 - Math.sqrt(det)) / (2.0 * a);
                }
            }
        }
    }

    private void fillLogS() { // Give a better name maybe? fillSurvProb?
        // We use the logS array to store the unlogged quantities first
        logS[0] = 1.0;

        // Solve integral equation with trapezoidal rule - corresponds to solving a system of linear equations
        for (int m = 1; m <= timeSteps; m++) {
            double trap_sum = 0.0;
            trap_sum += birthRateMat[0][m] * expmR[0][m] * p0[0] * logS[0] / 2.0;

            for (int i = 1; i < m; i++) {
                int j = m - i;
                trap_sum += birthRateMat[i][j] * expmR[i][j] * p0[i] * logS[i];
            }

            logS[m] = (2.0 * h * trap_sum + expmR[0][m]) / (1.0 - h * birthRateMat[m][0] * expmR[m][0] * p0[m]);
        }

        // Log and condition on non-extinction
        for (int i = 0; i <= timeSteps; i++) {
            logS[i] = Math.log(logS[i]) - Math.log(1.0 - p0[i]);
        }
    }

    private double branchLengthDensity(int k, int l) {
        double p0_k = p0[k];
        double one_minus_p0_k = 1.0 - p0_k;

        logBranchDens[0] = birthRateMat[k][0] * one_minus_p0_k * one_minus_p0_k;

        for (int m = 1; m <= l; m++) {
            double trap_sum = birthRateMat[k][m] * expmR[k][m] * p0_k * logBranchDens[0] / 2.0;

            for (int r = 1; r < m; r++) {
                int i = k + r;
                int j = m - r;

                trap_sum += birthRateMat[i][j] * expmR[i][j] * p0[i] * logBranchDens[r];
            }

            int i_end = k + m;
            double num = 2.0 * h * trap_sum + birthRateMat[k][m] * expmR[k][m] * one_minus_p0_k * one_minus_p0_k;
            double denom = 1.0 - h * birthRateMat[i_end][0] * expmR[i_end][0] * p0[i_end];

            logBranchDens[m] = num / denom;
        }

        return Math.log(logBranchDens[l]) - Math.log(1.0 - p0[k + l]);
    }

    private double calculateLogLikelihood() {
        logLikelihood = 0.0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);

            if (tree.isRoot(node)) {
                int k = (int) Math.round(tree.getNodeHeight(node) / h);
                int l = timeSteps - (int) Math.round(tree.getNodeHeight(node) / h);
                logLikelihood += branchLengthDensity(k, l);
            } else if (tree.isExternal(node)) {
                int l = (int) Math.round(tree.getBranchLength(node) / h);
                logLikelihood += logS[l];
            } else {
                int k = (int) Math.round(tree.getNodeHeight(node) / h);
                int l = (int) Math.round(tree.getBranchLength(node) / h);
                logLikelihood += branchLengthDensity(k, l);
            }
        }

        return logLikelihood;
    }



//    private double traversal(NodeRef node) { // Just some example code
//        if (tree.isExternal(node)) {
//            return tree.getBranchLength(node);
//        } else {
//            double leftCount = traversal(tree.getChild(node, 0));
//            double rightCount = traversal(tree.getChild(node, 1));
//            return leftCount + rightCount;
//        }
//    }
//
//    private double length() { // Just some example code
//        double length = 0;
//        for (int i = 0; i < tree.getNodeCount(); i++) {
//            length += tree.getBranchLength(tree.getNode(i));
//        }
//
//        return length;
//    }
}