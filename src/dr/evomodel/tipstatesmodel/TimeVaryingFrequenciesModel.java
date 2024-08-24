package dr.evomodel.tipstatesmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.TipStateAccessor;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

import static dr.util.Citation.Category.PRIOR_MODELS;

public class TimeVaryingFrequenciesModel extends AbstractModelLikelihood implements Citable {

    private final List<Epoch> epochs;

    @SuppressWarnings("unused")
    private final List<TipStateAccessor> accessors;
    private final Tree tree;
    private final Taxon taxon;

    @Override
    public Citation.Category getCategory() {
        return PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Time-varying integration of unobserved tip trait";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("MA", "Suchard"),
                    new Author("J", "Pekar"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static class Epoch {
        double cutOff;
        Parameter distribution;

        public Epoch(double cutOff, Parameter distribution) {
            this.cutOff = cutOff;
            this.distribution = distribution;
        }
    }

    public TimeVaryingFrequenciesModel(String name,
                                       List<TipStateAccessor> accessors,
                                       List<Epoch> epochs,
                                       Taxon taxon,
                                       Tree tree) {
        super(name);

        this.accessors = accessors;
        this.epochs = epochs;
        this.tree = tree;
        this.taxon = taxon;

        getTipNode(taxon);
    }

    public Tree getTree() { return tree; }

    public Taxon getTaxon() { return taxon; }

    public int getTipIndex(Taxon taxon) {
        int tipIndex = tree.getTaxonIndex(taxon);
        if (tipIndex == -1) {
            throw new RuntimeException("Unknown taxon");
        }
        return tipIndex;
    }

    private NodeRef getTipNode(Taxon taxon) {
        int tipIndex = getTipIndex(taxon);

        NodeRef tip = tree.getExternalNode(tipIndex);
        if (tree.getNodeTaxon(tip) != taxon) {
            throw new RuntimeException("Unknown taxon");
        }

        return tip;
    }

    private double calculateLogLikelihood() {
        return 0; // TODO Currently only valid for MH via TipStateOperator
    }

    public double[] getProbabilities(Taxon taxon) {

        double height = tree.getNodeHeight(getTipNode(taxon));

        int epochIndex = 0;
        while (height > epochs.get(epochIndex).cutOff) {
            ++epochIndex;
        }

        return epochs.get(epochIndex).distribution.getParameterValues();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new IllegalArgumentException("Unknown variable");
    }

    @Override
    public Model getModel() { return this; }

    @Override
    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }

    @Override
    public void makeDirty() { }
}
