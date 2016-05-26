/*
 * TreeTrait.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

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
     * @return Intent
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
     * Specifies whether this trait is loggable
     *
     * @return Intent
     */
    boolean getLoggable();

    /**
     * Default behavior
     */
    class DefaultBehavior {
        public boolean getLoggable() {
            return true;
        }

        public boolean getFormatAsArray() {
            return false;
        }
    }

    /**
     * An abstract base class for Double implementations
     */
    public abstract class D extends DefaultBehavior implements TreeTrait<Double> {

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
    public abstract class I extends DefaultBehavior implements TreeTrait<Integer> {

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
     * An abstract base class for String implementations
     */
    public abstract class S extends DefaultBehavior implements TreeTrait<String> {

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
    public abstract class DA extends DefaultBehavior implements TreeTrait<double[]> {

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
     * An abstract base class for String array implementations
     */
    public abstract class SA extends DefaultBehavior implements TreeTrait<String[]> {

        public Class getTraitClass() {
            return String[].class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return formatTrait(getTrait(tree, node), getFormatAsArray());
        }

        public static String formatTrait(String[] values, boolean asArray) {
            if (values == null || values.length == 0) return null;
            if (values.length > 1 || asArray) {
                StringBuilder sb = new StringBuilder("{");
                sb.append(values[0]);
                for (int i = 1; i < values.length; i++) {
                    sb.append(",");
                    sb.append(values[i]);
                }
                sb.append("}");

                return sb.toString();
            } else {
                return values[0];
            }
        }
    }

    /**
     * An abstract base class for int array implementations
     */
    public abstract class IA extends DefaultBehavior implements TreeTrait<int[]> {

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
    public abstract class SumOverTree<T> extends DefaultBehavior implements TreeTrait<T> {

        private static final String NAME_PREFIX = "sumOverTree_";
        private final TreeTrait<T> base;
        private final String name;
        private final boolean includeExternalNodes;
        private final boolean includeInternalNodes;

        public SumOverTree(TreeTrait<T> base) {
            this(NAME_PREFIX + base.getTraitName(), base);
        }

        public SumOverTree(String name, TreeTrait<T> base) {
            this(name, base, true, true);
        }

        public SumOverTree(String name, TreeTrait<T> base, boolean includeExternalNodes, boolean includeInternalNodes) {
            this.base = base;
            this.name = name;
            this.includeExternalNodes = includeExternalNodes;
            this.includeInternalNodes = includeInternalNodes;
        }

        public String getTraitName() {
            return name;
        }

        public Intent getIntent() {
            return Intent.WHOLE_TREE;
        }

        public T getTrait(Tree tree, NodeRef node) {
            T count = null;
            if (includeExternalNodes) {
                for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                    count = addToMatrix(count, base.getTrait(tree, tree.getExternalNode(i)));
                }
            }
            if (includeInternalNodes) {
                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                    count = addToMatrix(count, base.getTrait(tree, tree.getInternalNode(i)));
                }
            }
            return count;
        }

        public boolean getLoggable() {
            return base.getLoggable();
        }

        protected abstract T addToMatrix(T total, T summant);
    }

    /**
     * A wrapper class that sums a TreeTrait.DA over the entire tree
     */
    public class SumOverTreeDA extends SumOverTree<double[]> {

        public SumOverTreeDA(String name, TreeTrait<double[]> base,
                             boolean includeExternalNodes, boolean includeInternalNodes) {
            super(name, base, includeExternalNodes, includeInternalNodes);
        }

        public SumOverTreeDA(String name, TreeTrait<double[]> base) {
            super(name, base);
        }

        public SumOverTreeDA(TreeTrait<double[]> base) {
            super(base);
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return DA.formatTrait(getTrait(tree, node));
        }

        public Class getTraitClass() {
            return double[].class;
        }

        protected double[] addToMatrix(double[] total, double[] summant) {
            return addToMatrixStatic(total, summant);
        }

        protected static double[] addToMatrixStatic(double[] total, double[] summant) {
            if (summant == null) {
                return total;
            }
            final int length = summant.length;
            if (total == null) {
                total = new double[length];
            }
            for (int i = 0; i < length; i++) {
                total[i] += summant[i];
            }
            return total;
        }
    }

    /**
     * A wrapper class that sums a TreeTrait.D over the entire tree
     */
    public class SumOverTreeD extends SumOverTree<Double> {

        public SumOverTreeD(String name, TreeTrait<Double> base,
                            boolean includeExternalNodes, boolean includeInternalNodes) {
            super(name, base, includeExternalNodes, includeInternalNodes);
        }

        public SumOverTreeD(String name, TreeTrait<Double> base) {
            super(name, base);
        }

        public SumOverTreeD(TreeTrait<Double> base) {
            super(base);
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return D.formatTrait(getTrait(tree, node));
        }

        public Class getTraitClass() {
            return double[].class;
        }

        protected Double addToMatrix(Double total, Double summant) {
            return addToMatrixStatic(total, summant);
        }

        protected static Double addToMatrixStatic(Double total, Double summant) {
            if (summant == null) {
                return total;
            }
            if (total == null) {
                total = 0.0;
            }
            total += summant;
            return total;
        }
    }

    /**
     * An abstract wrapper class that sums a TreeTrait.Array into a TreeTrait
     */
    public abstract class SumAcrossArray<T, TA> extends DefaultBehavior implements TreeTrait<T> {

        private TreeTrait<TA> base;
        private String name;
        public static final String NAME_PREFIX = "sumAcrossArray_";

        public SumAcrossArray(TreeTrait<TA> base) {
            this(NAME_PREFIX + base.getTraitName(), base);
        }

        public SumAcrossArray(String name, TreeTrait<TA> base) {
            this.name = name;
            this.base = base;
        }

        public String getTraitName() {
            return name;
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

        public boolean getLoggable() {
            return base.getLoggable();
        }

        protected abstract T reduce(TA values);
    }

    /**
     * A wrapper class that sums a TreeTrait.DA into a TreeTrait.D
     */
    public class SumAcrossArrayD extends SumAcrossArray<Double, double[]> {

        public SumAcrossArrayD(String name, TreeTrait<double[]> base) {
            super(name, base);
        }

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


    /**
     * An abstract wrapper class that picks one entry out of TreeTrait<T> where T is an array
     */
    public abstract class PickEntry<T, TA> extends DefaultBehavior implements TreeTrait<T> {

        protected TreeTrait<TA> base;
        private String name;
        protected int index;

        public PickEntry(TreeTrait<TA> base, int index) {
            this(base.getTraitName() + "_" + (index + 1), base, index);
        }

        public PickEntry(String name, TreeTrait<TA> base, int index) {
            this.name = name;
            this.base = base;
//            if (base.getTraitClass() != int[].class || base.getTraitClass() != double[].class) {
//                throw new RuntimeException("Only supported for arrays");
//            }
            this.index = index;
        }

        public String getTraitName() {
            return name;
        }

        public Intent getIntent() {
            return base.getIntent();
        }
    }

    public class PickEntryD extends PickEntry<Double, double[]> {

        public PickEntryD(TreeTrait<double[]> base, int index) {
            super(base, index);
        }

        public PickEntryD(String name, TreeTrait<double[]> base, int index) {
            super(name, base, index);
        }

        public Class getTraitClass() {
            return Double.class;
        }

        public Double getTrait(Tree tree, NodeRef node) {
            return base.getTrait(tree, node)[index];
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return D.formatTrait(getTrait(tree, node));
        }
    }

    public class PickEntryDAndScale extends PickEntryD {

        public PickEntryDAndScale(TreeTrait<double[]> base, int index) {
            super(base, index);
        }

        public PickEntryDAndScale(String name, TreeTrait<double[]> base, int index) {
            super(name, base, index);
        }

        public Double getTrait(Tree tree, NodeRef node) {
            return (base.getTrait(tree, node)[index]) / tree.getBranchLength(node);
        }
    }

    public class PickEntryI extends PickEntry<Integer, int[]> {

        public PickEntryI(TreeTrait<int[]> base, int index) {
            super(base, index);
        }

        public PickEntryI(String name, TreeTrait<int[]> base, int index) {
            super(name, base, index);
        }

        public Class getTraitClass() {
            return Double.class;
        }

        public Integer getTrait(Tree tree, NodeRef node) {
            return base.getTrait(tree, node)[index];
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return I.formatTrait(getTrait(tree, node));
        }
    }

    /**
     * An abstract wrapper class that filters a TreeTrait<T>
     */
    public abstract class Filtered<T> extends DefaultBehavior implements TreeTrait<T> {

        private static final String NAME_PREFIX = "filtered_";
        private final TreeTrait<T> base;
        private final String name;
        private final TreeNodeFilter treeNodeFilter;

        public Filtered(TreeTrait<T> base, TreeNodeFilter treeNodeFilter) {
            this(NAME_PREFIX + base.getTraitName(), base, treeNodeFilter);
        }

        public Filtered(String name, TreeTrait<T> base, TreeNodeFilter treeNodeFilter) {
            this.base = base;
            this.name = name;
            this.treeNodeFilter = treeNodeFilter;
        }

        public String getTraitName() {
            return name;
        }

        public Intent getIntent() {
            return base.getIntent();
        }

        public T getTrait(Tree tree, NodeRef node) {

            T count = null;

            if (getIntent() == Intent.WHOLE_TREE) {
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    NodeRef tmpNode = tree.getNode(i);
                    if (treeNodeFilter.includeNode(tree, tmpNode)) {
                        count = addToMatrix(count, base.getTrait(tree, tmpNode));
                    }
                }
            } else {
                count = base.getTrait(tree, node);
                if (!treeNodeFilter.includeNode(tree, node)) {
                    count = zero(count);
                }
            }
            return count;
        }

        public boolean getLoggable() {
            return base.getLoggable();
        }

        protected abstract T addToMatrix(T total, T summant);

        protected abstract T zero(T copy);
    }

    public class FilteredD extends Filtered<Double> {

        public FilteredD(TreeTrait<Double> base, TreeNodeFilter treeNodeFilter) {
            super(base, treeNodeFilter);
        }

        public FilteredD(String name, TreeTrait<Double> base, TreeNodeFilter treeNodeFilter) {
            super(name, base, treeNodeFilter);
        }

        @Override
        protected Double addToMatrix(Double total, Double summant) {
            return SumOverTreeD.addToMatrixStatic(total, summant);
        }

        protected Double zero(Double copy) {
            return 0.0;
        }

        public Class getTraitClass() {
            return Double.class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return D.formatTrait(getTrait(tree, node));
        }
    }

    public class FilteredDA extends Filtered<double[]> {

        public FilteredDA(TreeTrait<double[]> base, TreeNodeFilter treeNodeFilter) {
            super(base, treeNodeFilter);
        }

        public FilteredDA(String name, TreeTrait<double[]> base, TreeNodeFilter treeNodeFilter) {
            super(name, base, treeNodeFilter);
        }

        @Override
        protected double[] addToMatrix(double[] total, double[] summant) {
            return SumOverTreeDA.addToMatrixStatic(total, summant);
        }

        protected double[] zero(double[] copy) {
            return new double[copy.length];
        }

        public Class getTraitClass() {
            return double[].class;
        }

        public String getTraitString(Tree tree, NodeRef node) {
            return DA.formatTrait(getTrait(tree, node));
        }
    }
}

