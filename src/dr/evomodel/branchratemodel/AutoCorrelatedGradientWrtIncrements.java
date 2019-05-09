package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class AutoCorrelatedGradientWrtIncrements implements GradientWrtParameterProvider, Reportable {

    private final AutoCorrelatedBranchRatesDistribution distribution;
    private final ArbitraryBranchRates branchRates;
    private final Tree tree;
    private Parameter parameter;
    private double[] cachedIncrements;

    public AutoCorrelatedGradientWrtIncrements(AutoCorrelatedBranchRatesDistribution distribution) {
        this.distribution = distribution;
        this.branchRates = distribution.getBranchRateModel();
        this.tree = distribution.getTree();
    }

    @Override
    public Likelihood getLikelihood() {
        return distribution.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        if (parameter == null) {
            parameter = createParameter();
        }
        return parameter;
    }

    @Override
    public int getDimension() {
        return distribution.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return distribution.getGradientWrtIncrements();
    }

    @Override
    public String getReport() {
        return new CheckGradientNumerically(
                this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null
        ).toString();
    }

    private Parameter createParameter() {
        return new Parameter.Proxy("increments", distribution.getDimension()) {

            @Override
            public double getParameterValue(int dim) {
                return distribution.getIncrement(dim);
            }

            @Override
            public void setParameterValue(int dim, double value) {
                throw new RuntimeException("Do not set single value at a time");
            }

            @Override
            public void setParameterValueQuietly(int dim, double value) {
                if (cachedIncrements == null) {
                    cachedIncrements = new double[getDimension()];
                }

                cachedIncrements[dim] = value;
            }

            @Override
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                throw new RuntimeException("Do not set single value at a time");
            }

            public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {

                double[] rates = new double[getDimension()];
                recurse(tree.getRoot(), rates, cachedIncrements, 1.0);

                Parameter rateParameter = distribution.getParameter();
                for (int i = 0; i < rates.length; ++i) {
                    rateParameter.setParameterValueQuietly(i, rates[i]);
                }

                rateParameter.fireParameterChangedEvent();
            }
        };
    }

    private void recurse(NodeRef node, double[] rates, double[] increments, double parentRate) {

        double rate = parentRate;

        if (!tree.isRoot(node)) {
            int index = branchRates.getParameterIndexFromNode(node);
            rate *= parentRate * distribution.inverseTransform(cachedIncrements[index]);

            rates[index] = rate;
        }

        if (!tree.isExternal(node)) {
            recurse(tree.getChild(node, 0), rates, increments, rate);
            recurse(tree.getChild(node, 1), rates, increments, rate);
        }
    }
}
