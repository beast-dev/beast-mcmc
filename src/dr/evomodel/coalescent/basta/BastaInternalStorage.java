package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.EigenDecomposition;

/**
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class BastaInternalStorage {

    double[] partials;
    double[] matrices;
    double[] coalescent;

    double[] storedMatrices;

    double[] e;
    double[] f;
    double[] g;
    double[] h;

    final double[] sizes;
    final double[] rates;

    final EigenDecomposition[] decompositions; // TODO flatten?

    private int currentNumCoalescentIntervals;
    private int currentNumPartials;

    final private int stateCount;

    public BastaInternalStorage(int maxNumCoalescentIntervals, int treeNodeCount, int stateCount) {

        this.currentNumPartials = 0;
        this.currentNumCoalescentIntervals = 0;

        this.stateCount = stateCount;

        this.sizes = new double[2 * stateCount];
        this.rates = new double[2 * stateCount];
        this.decompositions = new EigenDecomposition[1];

        resize(3 * treeNodeCount, maxNumCoalescentIntervals, null);
    }

    static private int getStartingPartialsCount(int maxNumCoalescentIntervals, int treeNodeCount) {
        return maxNumCoalescentIntervals * (treeNodeCount + 1); // TODO much too large
    }

    public void resize(int newNumPartials, int newNumCoalescentIntervals, BastaLikelihood likelihood) {

        if (newNumPartials > currentNumPartials) {
            this.partials = new double[newNumPartials * stateCount];
            this.currentNumPartials =  newNumPartials;
            if (likelihood != null) {
                likelihood.setTipData();
            }
        }

        if (newNumCoalescentIntervals > this.currentNumCoalescentIntervals) {
            this.matrices = new double[newNumCoalescentIntervals * stateCount * stateCount]; // TODO much too small (except for strict-clock)
            this.storedMatrices = new double[newNumCoalescentIntervals * stateCount * stateCount];
            this.coalescent = new double[newNumCoalescentIntervals];

            this.e = new double[newNumCoalescentIntervals * stateCount];
            this.f = new double[newNumCoalescentIntervals * stateCount];
            this.g = new double[newNumCoalescentIntervals * stateCount];
            this.h = new double[newNumCoalescentIntervals * stateCount];

            this.currentNumCoalescentIntervals = newNumCoalescentIntervals;
        }
    }

    public void storeState() {
        System.arraycopy(matrices, 0, storedMatrices, 0, currentNumCoalescentIntervals * stateCount * stateCount);
    }

    public void restoreState() {
        double[] temp = matrices;
        matrices = storedMatrices;
        storedMatrices = temp;
    }
}
