package dr.evomodel.tree;

import dr.inference.model.Parameter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Variable;

public class ExtendTipBranchTransform extends TreeTransform {
    public ExtendTipBranchTransform(Parameter diffusionPrecision) {
        super("extendTipBranchTransform");
        diffusionPrecision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.diffusionPrecision = diffusionPrecision;
        addVariable(diffusionPrecision);
    }

    @Override
    public double transform(Tree tree, NodeRef node, double originalHeight) {
        if (tree.isExternal(node)) {
            return originalHeight;
        }
        return originalHeight + getExtensionForNode(tree, node);
    }

    @Override
    protected double getScaleForNode(Tree tree, NodeRef node) {
        return diffusionPrecision.getParameterValue(0);
    }

    private double getExtensionForNode(Tree tree, NodeRef node) {
        double samplingVariance = 1.0 - 1.0 / getScaleForNode(tree, node);
        return samplingVariance * getScaleForNode(tree, node);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        fireModelChanged(diffusionPrecision);
    }

    @Override
    public String getInfo() {
        return null;
    }

    private final Parameter diffusionPrecision;
}
