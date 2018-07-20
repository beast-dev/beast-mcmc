package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface BranchSpecificFixedEffects {

    double getEffect(final Tree tree, final NodeRef node);

    ReadableVector getDesignVector(final Tree tree, final NodeRef node);

    Parameter getFixedEffectsParameter();

    class Default extends AbstractModel implements BranchSpecificFixedEffects {

        private final Parameter coefficients;
        private final List<CountableBranchCategoryProvider> categoryProviders;
        private final List<ContinuousBranchValueProvider> valueProviders;
        private final boolean hasIntercept = true;

        private final int dim;

        public Default(String name,
                       List<CountableBranchCategoryProvider> categoryProviders,
                       List<ContinuousBranchValueProvider> valueProviders,
                       Parameter coefficients) {
            super(name);

            this.coefficients = coefficients;
            this.categoryProviders = categoryProviders;
            this.valueProviders = valueProviders;

            this.dim = categoryProviders.size() + valueProviders.size() + (hasIntercept ? 1 : 0);

            if (coefficients.getDimension() != dim) {
                throw new IllegalArgumentException("Invalid parameter dimensions");
            }

            addModels(categoryProviders);
            addModels(valueProviders);
        }

        @Override
        public double getEffect(Tree tree, NodeRef node) {

            ReadableVector design = getDesignVector(tree, node);

            double sum = 0.0;
            for (int i = 0; i < design.getDim(); ++i) {
                sum += design.get(i) * coefficients.getParameterValue(i);
            }

            return sum;
        }

        @Override
        public ReadableVector getDesignVector(Tree tree, NodeRef node) {

            double[] design = new double[dim];

            int offset = 0;
            if (hasIntercept) {
                addIntercept(design);
                ++offset;
            }

            for (CountableBranchCategoryProvider categoryProvider : categoryProviders) {
                int category = categoryProvider.getBranchCategory(tree, node);
                design[category + offset] = 1.0;
            }
            offset += categoryProviders.size();

            for (ContinuousBranchValueProvider valueProvider : valueProviders) {
                design[offset] = valueProvider.getBranchValue(tree, node);
                ++offset;
            }

            return new WrappedVector.Raw(design);
        }

        private void addIntercept(double[] design) {
            design[0] = 1.0;
        }

        @Override
        public Parameter getFixedEffectsParameter() { return coefficients; }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) { }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }

        private void addModels(List list) {
            for (Object entry : list) {
                if (entry instanceof Model) {
                    addModel((Model) entry);
                }
            }
        }
    }
}
