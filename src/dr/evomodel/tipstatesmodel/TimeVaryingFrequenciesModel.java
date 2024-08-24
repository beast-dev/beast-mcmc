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
import java.util.Comparator;
import java.util.List;

import static dr.util.Citation.Category.PRIOR_MODELS;

public class TimeVaryingFrequenciesModel extends AbstractModelLikelihood implements Citable {

    private final List<Epoch> epochs;

    @SuppressWarnings("unused")
    private final List<TipStateAccessor> accessors;
    private final Tree tree;
    private final Taxon taxon;
    private final Parameter tipHeight;

    private boolean likelihoodKnown = false;
    private double storedLogLikelihood;
    private double logLikelihood;

    private static final boolean DEBUG = false;

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
                                       Tree tree,
                                       Parameter tipHeight) {
        super(name);

        this.accessors = accessors;
        this.epochs = epochs;
        this.tree = tree;
        this.taxon = taxon;
        this.tipHeight = tipHeight;

        if (tipHeight != null) {
            tipHeight.addVariableListener(this);
        }

        epochs.sort(Comparator.comparingDouble(lhs -> lhs.cutOff));

        getTipNode(taxon);
    }

    public List<TipStateAccessor> getAccessors() { return accessors; }

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

        if (DEBUG) {
            System.err.println("recompute prior");
        }

        int patternCount = accessors.get(0).getPatternCount();
        int[] currentStates = new int[patternCount];
        accessors.get(0).getTipStates(getTipIndex(taxon), currentStates);

        double[] probabilities = getProbabilities(taxon);
        int stateCount = probabilities.length;

        double logLikelihood = 0;
        for (int i = 0; i < patternCount; ++i) {
            int state = currentStates[i];
            if (state < stateCount) {
                logLikelihood += Math.log(probabilities[state]);
            }
        }

        return logLikelihood;
    }

    public double[] getProbabilities(Taxon taxon) {

        double height = tree.getNodeHeight(getTipNode(taxon));

//        if (DEBUG) {
//            System.err.println("height = " + height);
//        }

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
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
    }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == tipHeight) {
            likelihoodKnown = false;

            if (DEBUG) {
                System.err.println("height hit");
            }
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    public Model getModel() { return this; }

    @Override
    public double getLogLikelihood() {

        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        else if (DEBUG) {
            System.err.println("cached");
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;

//        if (DEBUG) {
//            System.err.println("make dirty");
//        }
    }
}
