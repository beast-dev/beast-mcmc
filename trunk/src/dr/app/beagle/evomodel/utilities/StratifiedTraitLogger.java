package dr.app.beagle.evomodel.utilities;

import dr.inference.loggers.Loggable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.evomodel.tree.TreeModel;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.StratifiedTraitOutputFormat;
import dr.evolution.tree.NodeRef;

/**
 * A utility class to format traits on a tree in a LogColumn-based file.  This class will take a TraitProvider
 * ultimately and its package location is bound to move
 *
 * @author Marc A. Suchard
 */

public class StratifiedTraitLogger implements Loggable {

    public StratifiedTraitLogger(String name,
                                 TreeModel tree,
                                 CodonPartitionedRobustCounting trait,
                                 StratifiedTraitOutputFormat format) {
        this.name = name;
        this.tree = tree;
        this.trait = trait;
        this.format = format;
    }

    public double[] getTotalTraitPerSite() {
        double[] count = new double[trait.getDimension()];

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                addToMatrix(count, trait.getExpectedCountsForBranch(node));  // TODO Change to function in TraitProvider
            }
        }
        return count;
    }

    public double getTotalTrait() {

        double[] count = getTotalTraitPerSite();
        double total = 0;
        for (double x : count) {
            total += x;
        }
        return total;
    }

    public double getTotalTraitForSite(int site) {

        if (site < 0 || site >= trait.getDimension()) {
            throw new RuntimeException("Invalid dimension #");
        }
        double total = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                total += trait.getExpectedCountsForBranch(node)[site]; // TODO Change to function in traitProvider
            }
        }
        return total;
    }

    private void addToMatrix(double[] total, double[] summant) {
        final int length = summant.length;
        for (int i = 0; i < length; i++) {
            total[i] += summant[i];
        }
    }

    public LogColumn[] getColumns() {

        if (format == StratifiedTraitOutputFormat.SUM_OVER_SITES) {
            return new LogColumn[]{
                    new CountColumn(name)
            };

        } else if (format == StratifiedTraitOutputFormat.PER_SITE ||
                format == StratifiedTraitOutputFormat.PER_SITE_WITH_UNCONDITIONED) {
            int nColumns = trait.getDimension();
            boolean allowUnconditioned = (format == StratifiedTraitOutputFormat.PER_SITE_WITH_UNCONDITIONED &&
                    (trait instanceof CodonPartitionedRobustCounting)); // TODO Replace with instance of UnconditionedTraitProvider

            if (allowUnconditioned) {
                nColumns++;
            }
            int index = 0;
            LogColumn[] allColumns = new LogColumn[nColumns];
            for (int j = 0; j < trait.getDimension(); j++) {
                allColumns[index++] = new ConditionedSiteCountColumn(name, j);
            }
            if (allowUnconditioned) {
                allColumns[index] = new UnconditionedSiteCountColumn(name);
            }
            return allColumns;

        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    private class CountColumn extends NumberColumn {

        public CountColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getTotalTrait();
        }
    }

    protected abstract class SiteCountColumn extends NumberColumn {

        protected int indexSite;

        public SiteCountColumn(String label, int j) {
            super(label + (j >= 0 ? "[" + (j + 1) + "]" : ""));
            indexSite = j;
        }
    }

    protected class ConditionedSiteCountColumn extends SiteCountColumn {

        public ConditionedSiteCountColumn(String label, int site) {
            super(label, site);
        }

        public double getDoubleValue() {
            return getTotalTraitForSite(indexSite);
        }
    }

    protected class UnconditionedSiteCountColumn extends SiteCountColumn {

        public UnconditionedSiteCountColumn(String label) {
            super("u_" + label, -1);
        }

        public double getDoubleValue() {
            return trait.getUnconditionedTraitValue(); // TODO Will check if instanceof UnconditionedTraitProvider
        }
    }

    private String name;
    private TreeModel tree;
    private CodonPartitionedRobustCounting trait; // TODO Change to TraitProvider
    private StratifiedTraitOutputFormat format;
}
