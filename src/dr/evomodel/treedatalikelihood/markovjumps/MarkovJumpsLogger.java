package dr.evomodel.treedatalikelihood.markovjumps;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.ProcessAlongTree;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;

public class MarkovJumpsLogger implements Loggable {

    public MarkovJumpsLogger(String tag,
                             ProcessAlongTree process,
                             TreeTrait trait,
                             TreeTrait marginalRateTrait,
                             TreeTrait unconditionalExpectedJumps,
                             TreeTrait siteSpecificRateTrait,
                             boolean reportUnconditional) {

        this.tag = tag;
        this.reportUnconditional = reportUnconditional;

        this.tree = process.getTree();

        Object value = trait.getTrait(tree, tree.getNode(0));
        if (value instanceof double[]) {
            this.trait = trait;
        } else {
            throw new IllegalArgumentException("Only double[] traits are currently supported");
        }

        this.traitDim = ((double[]) value).length;

        this.marginalRateTrait = marginalRateTrait;
        this.unconditionalExpectedJumps = unconditionalExpectedJumps;
        this.siteSpecificRateTrait = siteSpecificRateTrait;
    }

    @Override
    public LogColumn[] getColumns() {

        if (columns == null) {
            columns = new LogColumn[traitDim * (reportUnconditional ? 2 : 1)];

            for (int i = 0; i < traitDim; ++i) {
                columns[i] = new ConditionedCountColumn(tag, i);
            }

            if (reportUnconditional) {
                for (int i = 0; i < traitDim; ++i) {
                    columns[traitDim + i] = new UnconditionedCountColumn(tag, i);
                }
            }
        }

        return columns;
    }

    private abstract static class CountColumn extends NumberColumn {

        protected int indexSite;

        CountColumn(String label, int j) {
            super(label + (j >= 0 ? "[" + (j + 1) + "]" : ""));
            indexSite = j;
        }

        public abstract double getDoubleValue();
    }

    private class ConditionedCountColumn extends CountColumn {

        ConditionedCountColumn(String label, int j) {
            super("c_" + label, j);
        }

        public double getDoubleValue() {

            double total = 0;
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                if (tree.getRoot() != node) {
                    total += ((double[]) trait.getTrait(tree, node))[indexSite];
                }
            }

            return total;
        }
    }

    private class UnconditionedCountColumn extends CountColumn {

        public UnconditionedCountColumn(String label, int j) {
            super("u_" + label, j);
        }

        public double getDoubleValue() {



            double[] siteRates = (double[]) siteSpecificRateTrait.getTrait(tree, null);
            double[] marginalRates = (double[]) marginalRateTrait.getTrait(tree, null);

            double total = 0;
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                if (tree.getRoot() != node) {
                    final int substitutionModelIndex = 0;
                    double expectedLength = ((double[]) unconditionalExpectedJumps.getTrait(tree, node))[indexSite];
                    // was:  double value = markovjumps.get(indexRegistration).getMarginalRate() * getExpectedTreeLength() * siteRateModel.getRateForCategory(rateCategory[indexSite])
                    total += marginalRates[substitutionModelIndex] * expectedLength;
                }
            }
            total *= siteRates[indexSite];

            return total;
        }

    }

    private LogColumn[] columns;

    private final String tag;
    private final Tree tree;
    private final TreeTrait trait;
    private final TreeTrait marginalRateTrait;
    private final TreeTrait unconditionalExpectedJumps;
    private final TreeTrait siteSpecificRateTrait;
    private final int traitDim;
    private final boolean reportUnconditional;
}
