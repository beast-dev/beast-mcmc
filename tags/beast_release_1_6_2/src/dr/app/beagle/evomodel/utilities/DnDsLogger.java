package dr.app.beagle.evomodel.utilities;

import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.EmpiricalBayesPoissonSmoother;

/**
 * @author Philippe Lemey
 * @author Marc A. Suchard
 */
public class DnDsLogger implements Loggable {

    public DnDsLogger(String name, Tree tree, TreeTrait[] traits, boolean useSmoothing) {
        this.tree = tree;
        this.traits = traits;
        numberSites = getNumberSites();
        this.name = name;
        this.useSmoothing = useSmoothing;

        for (int i = 0; i < NUM_TRAITS; i++) {
            if (traits[i].getIntent() != TreeTrait.Intent.WHOLE_TREE) {
                System.err.println("BAD");
                System.exit(-1);
            }
        }
    }

    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[numberSites];
        for (int i = 0; i < numberSites; i++) {
            columns[i] = new SmoothedColumn(name, i);
        }
        return columns;
    }

    private class SmoothedColumn extends NumberColumn {

        private final int index;

        public SmoothedColumn(String label, int index) {
            super(label + "[" + (index+1) + "]");
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            if (index == 0) { // Assumes that columns are accessed IN ORDER
                doSmoothing();
            }
            return doCalculation(index);
        }
    }

    private double doCalculation(int index) {
        return (cachedValues[CN][index] / cachedValues[UN][index]) /
                (cachedValues[CS][index] / cachedValues[US][index]);
    }

    private int getNumberSites() {
        double[] values = (double[]) traits[0].getTrait(tree, tree.getRoot());
        return values.length;
    }

    private void doSmoothing() {

        if (cachedValues == null) {
            cachedValues = new double[NUM_TRAITS][];
        }

        for (int i = 0; i < NUM_TRAITS; i++) {
            cachedValues[i] = (double[]) traits[i].getTrait(tree, tree.getRoot());
            if (useSmoothing) {
                cachedValues[i] = EmpiricalBayesPoissonSmoother.smooth((double[]) traits[i].getTrait(tree, tree.getRoot()));
            }
        }
    }

    private final TreeTrait[] traits;
    private final Tree tree;
    private final int numberSites;
    private final String name;
    private final boolean useSmoothing;

    private final static int NUM_TRAITS = 4;
    private final static int CS = 0;
    private final static int US = 1;
    private final static int CN = 2;
    private final static int UN = 3;

    private double[][] cachedValues;

    public static String[] traitNames = new String[] {
            CodonPartitionedRobustCounting.SITE_SPECIFIC_PREFIX + CodonLabeling.SYN.getText(),
            CodonPartitionedRobustCounting.UNCONDITIONED_PREFIX + CodonLabeling.SYN.getText(),
            CodonPartitionedRobustCounting.SITE_SPECIFIC_PREFIX + CodonLabeling.NON_SYN.getText(),
            CodonPartitionedRobustCounting.UNCONDITIONED_PREFIX + CodonLabeling.NON_SYN.getText()
    };
}


