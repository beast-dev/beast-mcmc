package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.MathUtils;

import java.util.Arrays;

public class NewBDSSHistorySimulator implements Loggable {

    public NewBDSSHistorySimulator(NewBirthDeathSerialSamplingModel bdss,
                                   Tree tree,
                                   boolean startAtOrigin,
                                   boolean conditionOnSurvivalToPresent,
                                   boolean recordBeforeSampling) {
        this.bdss = bdss;
        this.tree = tree;
        this.startAtOrigin = startAtOrigin;
        this.conditionOnSurvivalToPresent = conditionOnSurvivalToPresent;
        this.recordBeforeSampling = recordBeforeSampling;
        this.dim = bdss.getDeathRateParameter().getDimension();

        history = new double[dim];
        nSeen = 0;
    }

    private double getStartTime() {
        double t = bdss.originTime.getParameterValue(0);
        if ( !startAtOrigin ) {
            t = tree.getNodeHeight(tree.getRoot());
        }
        return t;
    }

    private int getStartingIndex() {
        double t = getStartTime();
        double[] grid = bdss.getBreakPoints();
        int i = 0;
        while (grid[i] < t) {
            i++;
        }
        return i;
    }

    private int drawNext(int idx, int n0, double tOld, double tYoung) {
        double birth = bdss.getBirthRateParameter().getParameterValue(idx);
        double death = bdss.getDeathRateParameter().getParameterValue(idx);
        double sampling = bdss.getSamplingRateParameter().getParameterValue(idx);
        double treatment = bdss.getTreatmentProbabilityParameter().getParameterValue(idx);

        double totalDeath = death + sampling * treatment;
        double netBirth = birth - totalDeath;

        int n = 0;

        double expRT = Math.exp(netBirth * (tOld - tYoung));

        // We are assuming birthRate != totalDeathRate
        double surivalProb = 1.0 - (totalDeath * (expRT - 1.0))/(birth * expRT - totalDeath);
        double beta = birth/totalDeath * (1.0 - surivalProb);
        double rate = -Math.log(beta);
        for (int i = 0; i < n0; i++) {
            // Tavare 2018, doi:10.1017/apr.2018.84, Eqns 10 and 4
            if (MathUtils.nextDouble() < surivalProb) {
                n += 1 + Math.floor(MathUtils.nextExponential(rate));
            }
        }

        return n;
    }

    private int doMassSampling(int idx, int n0) {
        int n = n0;
        double rho = bdss.getSamplingProbabilityParameter().getParameterValue(idx);
        double treatment = bdss.getTreatmentProbabilityParameter().getParameterValue(idx);

        if (rho > 0.0) {
            for (int i = 0; i < n; i++) {
                if (MathUtils.nextDouble() < rho * treatment) {
                    n--;
                }
            }
        }

        return n;
    }

    private void simulateHistory() {
        boolean done = false;

        while (!done) {
            Arrays.fill(history,0.0);

            double t = getStartTime();
            int index = getStartingIndex();
            int n = startAtOrigin ? 1 : 2;

            double[] grid = bdss.getBreakPoints();

            while (index > -1) {
                double gridIntervalStart = (index == 0) ? 0.0 : grid[index - 1];
                n = drawNext(index, n, t, gridIntervalStart);
                history[index] = (double)n;
                n = doMassSampling(index,n);
                if (!recordBeforeSampling) {
                    history[index] = (double)n;
                }
                t = gridIntervalStart;
                index--;
            }
            if (conditionOnSurvivalToPresent == false || history[0] > 0) {
                done = true;
            }
        }
    }

    private double getHistory(int i) {
        return history[i];
    }

    private class NewBDSSHistorySimulatorColumn extends NumberColumn {

        final int index;

        public NewBDSSHistorySimulatorColumn(int index) {
            super("NewBDSSHistorySimulatorColumn");
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            // Poor man's dirty/clean
            if (index == 0) {
                simulateHistory();
            }
            return getHistory(index);
        }

    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = new NewBDSSHistorySimulatorColumn[dim];
            for (int index = 0; index < dim; index++) {
                columns[index] = new NewBDSSHistorySimulatorColumn(index);
            }
        }

        return columns;
    }

    private final NewBirthDeathSerialSamplingModel bdss;
    private final Tree tree;
    private final boolean startAtOrigin;
    private final boolean conditionOnSurvivalToPresent;
    private final boolean recordBeforeSampling;
    private final int dim;
    private double[] history;
    private int nSeen;
    private NewBDSSHistorySimulatorColumn[] columns = null;

}