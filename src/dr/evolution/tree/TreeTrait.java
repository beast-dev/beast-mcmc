package dr.evolution.tree;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @author Alexei Drummond
 * @version $Id$
 */
public interface TreeTrait<T> {
    public enum Intent {
        NODE,
        BRANCH,
        WHOLE_TREE
    }

    /**
     * The human readable name of this trait
     *
     * @return the name
     */
    String getTraitName();

    /**
     * Specifies whether this is a trait of the tree, the nodes or the branch
     *
     * @return true if a branch property
     */
    Intent getIntent();

    /**
     * Return a class object for the trait
     *
     * @return the class
     */
    Class getTraitClass();

    /**
     * Returns the trait values for the given node. If this is a branch trait then
     * it will be for the branch above the specified node (and may not be valid for
     * the root). The array will be the length returned by getDimension().
     *
     * @param tree a reference to a tree
     * @param node a reference to a node
     * @return the trait value
     */
    T getTrait(final Tree tree, final NodeRef node);

    /**
     * Get a string representations of the trait value.
     *
     * @param tree a reference to a tree
     * @param node a reference to a node
     * @return the trait string representation
     */
    String getTraitString(final Tree tree, final NodeRef node);

    /**
     * An abstract base class for Double implementations
     */
    public abstract class D implements TreeTrait<Double> {

        public Class getTraitClass() {
            return Double.class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return formatTrait(getTrait(tree, node));
        }

        public static String formatTrait(Double value) {
            if (value == null) {
                return null;
            }
            return value.toString();
        }
    }

    /**
     * An abstract base class for Double implementations
     */
    public abstract class I implements TreeTrait<Integer> {

        public Class getTraitClass() {
            return Integer.class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return formatTrait(getTrait(tree, node));
        }

        public static String formatTrait(Integer value) {
            if (value == null) {
                return null;
            }
            return value.toString();
        }
    }

    /**
     * An abstract base class for Double implementations
     */
    public abstract class S implements TreeTrait<String> {

        public Class getTraitClass() {
            return String.class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return getTrait(tree, node);
        }
    }

    /**
     * An abstract base class for double array implementations
     */
    public abstract class DA implements TreeTrait<double[]> {

        public Class getTraitClass() {
            return double[].class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return formatTrait(getTrait(tree, node));
        }

        public static String formatTrait(double[] values) {
            if (values == null || values.length == 0) return null;
            if (values.length > 1) {
                StringBuilder sb = new StringBuilder("{");
                sb.append(values[0]);
                for (int i = 1; i < values.length; i++) {
                    sb.append(",");
                    sb.append(values[i]);
                }
                sb.append("}");

                return sb.toString();
            } else {
                return Double.toString(values[0]);
            }
        }
    }

    /**
     * An abstract base class for int array implementations
     */
    public abstract class IA implements TreeTrait<int[]> {

        public Class getTraitClass() {
            return int[].class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return formatTrait(getTrait(tree, node));
        }

        public static String formatTrait(int[] values) {
            if (values == null || values.length == 0) return null;
            if (values.length > 1) {
                StringBuilder sb = new StringBuilder("{");
                sb.append(values[0]);
                for (int i = 1; i < values.length; i++) {
                    sb.append(",");
                    sb.append(values[i]);
                }
                sb.append("}");

                return sb.toString();
            } else {
                return Integer.toString(values[0]);
            }
        }
    }

    /**
     * An abstract wrapper class that sums a TreeTrait<T> over the entire tree
     */
    public abstract class SumOverTree<T> implements TreeTrait<T> {

        private static final String NAME_PREFIX = "";
        private TreeTrait<T> base;

        public SumOverTree(TreeTrait<T> base) {
            this.base = base;
        }

        public String getTraitName() {
            return NAME_PREFIX + base.getTraitName();
        }

        public Intent getIntent() {
            return Intent.WHOLE_TREE;
        }

        public T getTrait(Tree tree, NodeRef node) {
            T count = null;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                addToMatrix(count, base.getTrait(tree, tree.getNode(i)));
            }
            return count;
        }

        protected abstract void addToMatrix(T total, T summant);
    }

    /**
     * A wrapper class that sums a TreeTrait.DA over the entire tree
     */
    public class SumOverTreeDA extends SumOverTree<double[]> {

        public SumOverTreeDA(TreeTrait<double[]> base) {
            super(base);
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return DA.formatTrait(getTrait(tree, node));
        }

        public Class getTraitClass() {
            return double[].class;
        }

        protected void addToMatrix(double[] total, double[] summant) {
            if (summant == null) {
                return;
            }
            final int length = summant.length;
            if (total == null) {
                total = new double[length];
            }
            for (int i = 0; i < length; i++) {
                total[i] += summant[i];
            }
        }
    }

    /**
     * A wrapper class that sums a TreeTrait.D over the entire tree
     */
    public class SumOverTreeD extends SumOverTree<Double> {

        public SumOverTreeD(TreeTrait<Double> base) {
            super(base);
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return D.formatTrait(getTrait(tree, node));
        }

        public Class getTraitClass() {
            return int[].class;
        }

        protected void addToMatrix(Double total, Double summant) {
            if (summant == null) {
                return;
            }
            if (total == null) {
                total = 0.0;
            }
            total += summant;
        }
    }

    /**
     * An abstract wrapper class that sums a TreeTrait.Array into a TreeTrait
     */
    public abstract class SumAcrossArray<T,TA> implements TreeTrait<T> {

        private TreeTrait<TA> base;
        public static final String NAME_PREFIX = "";

        public SumAcrossArray(TreeTrait<TA> base) {
            this.base = base;
        }

        public String getTraitName() {
            return NAME_PREFIX + base.getTraitName();
        }

        public Intent getIntent() {
            return base.getIntent();
        }

        public T getTrait(Tree tree, NodeRef node) {
            TA values = base.getTrait(tree, node);
            if (values == null) {
                return null;
            }
            return reduce(values);

        }

        protected abstract T reduce(TA values);
    }

    /**
     * A wrapper class that sums a TreeTrait.DA into a TreeTrait.D
     */
    public class SumAcrossArrayD extends SumAcrossArray<Double, double[]> {

        public SumAcrossArrayD(TreeTrait<double[]> base) {
            super(base);
        }

        public Class getTraitClass() {
            return Double.class;
        }

        protected Double reduce(double[] values) {
            double total = 0.0;
            for (double value : values) {
                total += value;
            }
            return total;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return D.formatTrait(getTrait(tree, node));
        }
    }
}

