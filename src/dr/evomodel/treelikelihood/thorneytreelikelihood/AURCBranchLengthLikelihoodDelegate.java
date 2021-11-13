package dr.evomodel.treelikelihood.thorneytreelikelihood;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.inference.model.*;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.NegativeBinomialDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

// TODO check if this can be used as a true branch rate model?
public class AURCBranchLengthLikelihoodDelegate extends AbstractBranchRateModel implements ThorneyBranchLengthLikelihoodDelegate, Citable {
    private final Parameter mu;
    private final Parameter omega;
    private final double scale; // for rate statistic

    public AURCBranchLengthLikelihoodDelegate(String name, Parameter mu, Parameter omega, double scale) {
        super(name);
        this.mu = mu;
        addVariable(mu);
        this.omega = omega;
        addVariable(omega);
        this.scale = scale;
    }

    @Override
    public double getLogLikelihood(double observed, Tree tree, NodeRef node) {
        double time = tree.getBranchLength(node);

        return NegativeBinomialDistribution.logPdf(observed,
                this.mu.getParameterValue(0) * time / this.omega.getParameterValue(0),
                this.omega.getParameterValue(0) / (1 + this.omega.getParameterValue(0)));
    }

    @Override
    public double getGradientWrtTime(double mutations, double time) {
        throw new RuntimeException("Not implemented!");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        // just samples from gamma distribution
        double mu = this.mu.getParameterValue(0);
        double omega = this.omega.getParameterValue(0);
        double l = tree.getBranchLength(node);
        double shape = mu*l/omega;
        double scale =omega/l;

        return GammaDistribution.nextGamma(shape, scale)/this.scale;

    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        return "Additive Uncorrelated Relaxed Clock model";
    }

    /**
     * @return a list of citations associated with this object
     */
    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }
    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Didelot"),
                    new Author("I", "Siveroni"),
                    new Author("EM", "Volz")
            },
            "Additive uncorrelated relaxed clock models for the dating of genomic epidemiology phylogenies",
            2020,
            "Mol Biol Evol",
            38,307,317

    );
}
