package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.ejml.data.DenseMatrix64F;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
abstract class TransformedTreeTrait extends AbstractModel implements TreeTrait<double[]> {

    private TransformedTreeTrait(TreeTrait<double[]> baseTrait, List<Tree> trees) {
        super(getName(baseTrait.getTraitName()));

        for (Tree tree : trees) {
            if (tree instanceof TreeModel) {
                addModel((TreeModel) tree);
            }
        }

        this.baseTrait = baseTrait;
        this.trees = trees;
        this.transformKnown = false;
    }

    @Override
    public String getTraitName() {
        return getName(baseTrait.getTraitName());
    }

    @Override
    public Intent getIntent() {
        return baseTrait.getIntent();
    }

    @Override
    public Class getTraitClass() {
        return baseTrait.getTraitClass();
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {

        assert (trees.contains(tree));

        if (!transformKnown) {
            updateTransform();
            transformKnown = true;
        }

        return transform(baseTrait.getTrait(tree, node));
    }

    @Override
    public String getTraitString(Tree tree, NodeRef node) {

        assert (trees.contains(tree));

        return TreeTrait.DA.formatTrait(getTrait(tree, node));
    }

    @Override
    public boolean getLoggable() {
        return baseTrait.getLoggable();
    }

    public class MultiDimensionalScaling extends TransformedTreeTrait {

        public MultiDimensionalScaling(TreeTrait<double[]> baseTrait, List<Tree> trees) {
            this(baseTrait, trees, trees.get(0));
        }

        public MultiDimensionalScaling(TreeTrait<double[]> baseTrait, List<Tree> trees,
                                       Tree sentinelTree) {
            super(baseTrait, trees);

            this.sentinelTree = sentinelTree;
            this.dim = baseTrait.getTrait(sentinelTree, sentinelTree.getExternalNode(0)).length;
        }

        void updateTransform() {
            DenseMatrix64F sentinel = getSentinelValue();

            linearOperator = null; // TODO use firstValue
        }

        private DenseMatrix64F getSentinelValue() {

            double[][] value = new double[dim][];

            for (int i = 0; i < dim; ++i) {
                value[i] = baseTrait.getTrait(sentinelTree, sentinelTree.getExternalNode(i));
            }

            return new DenseMatrix64F(value);
        }

        private final Tree sentinelTree;
        private final int dim;
    }

    public class PhylogeneticFactorAnalysis extends TransformedTreeTrait {

        public PhylogeneticFactorAnalysis(TreeTrait<double[]> baseTrait, List<Tree> trees, Parameter loadings) {
            super(baseTrait, trees);

            addVariable(loadings);

            this.loadings = loadings;
        }

        void updateTransform() {
              linearOperator = null; // TODO use firstValue
        }
        
        private final Parameter loadings;
    }
    
    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        transformKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        transformKnown = false;
    }

    private double[] transform(double[] value) {

        final int dim = value.length;

        double[] result = new double[dim];

        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;

            for (int j = 0; j < dim; ++j) {
                sum += linearOperator.unsafe_get(i, j);
            }

            result[i] = sum;
        }

        return result;
    }

    abstract void updateTransform();

    private static String getName(String name) {
        return "transformed." + name;
    }

    private final TreeTrait<double[]> baseTrait;
    private final List<Tree> trees;
    private boolean transformKnown;

    private DenseMatrix64F linearOperator;
}
