package dr.evomodel.branchratemodel;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.Transform;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface BranchSpecificFixedEffects {

    double getEffect(final Tree tree, final NodeRef node);

    double[] getDesignVector(final Tree tree, final NodeRef node);

    Parameter getFixedEffectsParameter();

    double[] getDifferential(double rate, final Tree tree, final NodeRef node);

    int getDimension();

    abstract class Base extends AbstractModel implements BranchSpecificFixedEffects {

        public Base(String name) {
            super(name);
        }

        public double[] getDifferential(double rate, final Tree tree, final NodeRef node) {
            double[] result = getDesignVector(tree, node);
            final double multiplier = rate / getEffect(tree, node);
            for (int i = 0; i < result.length; i++) {
                result[i] *= multiplier;
            }
            return result;
        }
    }

    class None extends Base implements BranchSpecificFixedEffects {

        private final Parameter location;
        private final static double[] design = new double[] { 1.0 };

        public None(Parameter location) {
            super("No effects");
            this.location = location;

            addVariable(location);
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            // Do nothing
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == location) {
                fireModelChanged();
            } else {
                throw new RuntimeException("Unknown variable: " + variable.getVariableName());
            }
        }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }

        @Override
        public double getEffect(Tree tree, NodeRef node) {
            return location.getParameterValue(0);
        }

        @Override
        public double[] getDesignVector(Tree tree, NodeRef node) {
            return design.clone();
        }

        @Override
        public Parameter getFixedEffectsParameter() {
            return location;
        }

        @Override
        public int getDimension() {
            return 1;
        }
    }

    class Transformed extends Base implements BranchSpecificFixedEffects {

        private final BranchSpecificFixedEffects effects;
        private final Transform transform;

        public Transformed(BranchSpecificFixedEffects effects, Transform transform) {
            super("With transform");
            this.effects = effects;
            this.transform = transform;
            addModel((Model) effects);
        }

        @Override
        public double[] getDifferential(double rate, final Tree tree, final NodeRef node) {
            double[] result = super.getDifferential(rate, tree, node);
            final double multiplier = transform.gradient(getEffect(tree, node));
            for (int i = 0; i < result.length; i++) {
                result[i] *= multiplier;
            }
            return result;
         }

        @Override
        public int getDimension() {
            return effects.getDimension();
        }

        @Override
        public double getEffect(Tree tree, NodeRef node) {
            double transformedEffect = effects.getEffect(tree, node);
            return transform.inverse(transformedEffect);
        }

        @Override
        public double[] getDesignVector(Tree tree, NodeRef node) {
            return effects.getDesignVector(tree, node);
        }

        @Override
        public Parameter getFixedEffectsParameter() {
            return effects.getFixedEffectsParameter();
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == effects.getFixedEffectsParameter()) {
                fireModelChanged();
            } else {
                throw new RuntimeException("Unknown variable: " + variable.getVariableName());
            }
        }

        @Override
        protected void storeState() {

        }

        @Override
        protected void restoreState() {

        }

        @Override
        protected void acceptState() {

        }
    }

    class Default extends Base implements BranchSpecificFixedEffects, Citable {

        private final Parameter coefficients;
        private final List<CountableBranchCategoryProvider> categoryProviders;
        private final List<ContinuousBranchValueProvider> valueProviders;
        private final List<BranchRates> branchRateProviders;
        private final boolean hasIntercept;

        private final int dim;

        public Default(String name,
                       List<CountableBranchCategoryProvider> categoryProviders,
                       List<ContinuousBranchValueProvider> valueProviders,
                       List<BranchRates> branchRateProviders,
                       Parameter coefficients,
                       boolean hasIntercept) {
            super(name);

            this.coefficients = coefficients;
            this.categoryProviders = categoryProviders;
            this.valueProviders = valueProviders;
            this.branchRateProviders = branchRateProviders;
            this.hasIntercept = hasIntercept;

            this.dim = categoryProviders.size() +
                    valueProviders.size() +
                    branchRateProviders.size() +
                    (hasIntercept ? 1 : 0);  

            if (coefficients.getDimension() != dim) {
                throw new IllegalArgumentException("Invalid parameter dimensions");
            }

            addModels(categoryProviders);
            addModels(valueProviders);
            addModels(branchRateProviders);
            addVariable(coefficients);
        }

        @Override
        public double getEffect(Tree tree, NodeRef node) {

            double[] design = getDesignVector(tree, node);

            double sum = 0.0;
            for (int i = 0; i < dim; ++i) {
                sum += design[i] * coefficients.getParameterValue(i);
            }

            return sum;
        }

        @Override
        public double[] getDesignVector(Tree tree, NodeRef node) {

            double[] design = new double[dim];

            int offset = 0;
            if (hasIntercept) {
                addIntercept(design);
                ++offset;
            }

            for (CountableBranchCategoryProvider categoryProvider : categoryProviders) {
                int category = categoryProvider.getBranchCategory(tree, node);
                if (category != 0) {
                    design[(category - 1) + offset] = 1.0;
                }
            }
            offset += categoryProviders.size();

            for (ContinuousBranchValueProvider valueProvider : valueProviders) {
                design[offset] = valueProvider.getBranchValue(tree, node);
                ++offset;
            }

            for (BranchRates branchRates : branchRateProviders) {
                design[offset] = transformFromBranchRateModel(branchRates.getBranchRate(tree, node));
                ++offset;
            }

            return design;
        }

        private double transformFromBranchRateModel(double x) {
            return Math.log(x);
        }

        private void addIntercept(double[] design) {
            design[0] = 1.0;
        }

        @Override
        public Parameter getFixedEffectsParameter() { return coefficients; }

        @Override
        public int getDimension() {
            return dim;
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == getFixedEffectsParameter()) {
                fireModelChanged();
            } else {
                throw new RuntimeException("Unknown variable: " + variable.getVariableName());
            }
        }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }

        public double[][] getDesignMatrix(Tree tree) {

            double[][] matrix = new double[tree.getNodeCount() - 1][];

            int offset = 0;
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                if (node != tree.getRoot()) {
                    matrix[offset] = getDesignVector(tree, node);
                    ++offset;
                }
            }
            return matrix;
        }

        private void addModels(List list) {
            for (Object entry : list) {
                if (entry instanceof Model) {
                    addModel((Model) entry);
                }
            }
        }

        @Override
        public Citation.Category getCategory() {
            return Citation.Category.MOLECULAR_CLOCK;
        }

        @Override
        public String getDescription() {
            return "Location-scale relaxed clock";
        }

        @Override
        public List<Citation> getCitations() {
            return Collections.singletonList(CITATION);
        }

        public static Citation CITATION = new Citation(
                new Author[]{
                        new Author("X", "Ji"),
                        new Author("P", "Lemey"),
                        new Author("MA", "Suchard")
                },
                Citation.Status.IN_PREPARATION
        );
    }
}
