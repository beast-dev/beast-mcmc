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
     * The number of dimensions for this trait
     * @return the trait dimension
     */
    int getDimension();

    /**
     * Returns an array giving the trait values for the given node. If this
     * is a branch trait then it will be for the branch above the specified
     * node (and may not be valid for the root). The array will be the length
     * returned by getDimension().
     * @param tree
     * @param node a reference to a node
     * @return the trait values
     */
    T[] getTrait(final Tree tree, final NodeRef node);

    /**
     * Get an array containing string representations of the trait values.
     * @param tree
     * @param node a reference to a node
     * @return the trait string representations
     */
    String[] getTraitString(final Tree tree, final NodeRef node);

     /**
     * An abstract base class for Double implementations
     */
    public abstract class D implements TreeTrait<Double> {

        public Class getTraitClass() {
            return Double.class;
        }

        public String[] getTraitString(Tree tree, NodeRef node) {
            Double[] values = getTrait(tree, node);
            if (values == null) {
                return null;
            }
            String[] strings = new String[values.length];
            for (int i = 0; i < values.length; i++) {
                strings[i] = values[i].toString();
            }
            return strings;
        }

         /**
          * A static utility function to convert an array of native doubles into
          * an equivalent array of Double objects. 
          * @param values
          * @return
          */
        public static Double[] toArray(double[] values) {
            Double[] array = new Double[values.length];
            for (int i = 0; i < values.length; i++) {
                array[i] = values[i];
            }
            return array;
        }
    }

    /**
     * An abstract base class for Double implementations
     */
    public abstract class I implements TreeTrait<Integer> {

        public Class getTraitClass() {
            return Double.class;
        }

        public String[] getTraitString(Tree tree, NodeRef node) {
            Integer[] values = getTrait(tree, node);
            if (values == null) {
                return null;
            }
            String[] strings = new String[values.length];
            for (int i = 0; i < values.length; i++) {
                strings[i] = values[i].toString();
            }
            return strings;
        }

        /**
         * A static utility function to convert an array of native ints into
         * an equivalent array of Integer objects.
         * @param values
         * @return
         */
        public static Integer[] toArray(int[] values) {
            Integer[] array = new Integer[values.length];
            for (int i = 0; i < values.length; i++) {
                array[i] = values[i];
            }
            return array;
        }
    }

    /**
     * An abstract base class for Double implementations
     */
    public abstract class S implements TreeTrait<String> {

        public Class getTraitClass() {
            return String.class;
        }

        public String[] getTraitString(Tree tree, NodeRef node) {
            return getTrait(tree, node);
        }
    }
}

