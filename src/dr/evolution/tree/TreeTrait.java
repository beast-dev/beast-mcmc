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
        BRANCH
    }

    /**
     * The human readable name of this trait
     * @return the name
     */
    String getTraitName();

    /**
     * Specifies whether this is a trait of the tree, the nodes or the branch
     * @return true if a branch property
     */
    Intent getIntent();

    /**
     * Return a class object for the trait
     * @return the class
     */
    Class getTraitClass();

    /**
     * Returns the trait values for the given node. If this is a branch trait then
     * it will be for the branch above the specified node (and may not be valid for
     * the root). The array will be the length returned by getDimension().
     * @param tree a reference to a tree
     * @param node a reference to a node
     * @return the trait value
     */
    T getTrait(final Tree tree, final NodeRef node);

    /**
     * Get a string representations of the trait value.
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
            Double value = getTrait(tree, node);
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
            return Double.class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            Integer value = getTrait(tree, node);
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
            double[] values = getTrait(tree, node);
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
            int[] values = getTrait(tree, node);
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

}

