/*
 * MascotLikelihood.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

/**
 * BEAST-X model wrapper for {@link MascotCore}.
 */
public class MascotLikelihood extends AbstractModelLikelihood implements Reportable {

    public static final String MASCOT_LIKELIHOOD = "mascotLikelihood";

    private final TreeModel treeModel;
    private final Parameter tipStates;
    private final MascotDynamics dynamics;
    private final double maxStep;
    private final boolean checkProbabilities;

    private MascotEventTape eventTape;
    private MascotCore core;
    private boolean eventTapeKnown;
    private boolean coreKnown;
    private boolean likelihoodKnown;

    private double logLikelihood;
    private double storedLogLikelihood;

    public MascotLikelihood(String name,
                            TreeModel treeModel,
                            Parameter tipStates,
                            Parameter theta,
                            Parameter epochTimes,
                            int stateCount,
                            double maxStep,
                            boolean checkProbabilities) {
        super(name == null ? MASCOT_LIKELIHOOD : name);
        this.treeModel = treeModel;
        this.tipStates = tipStates;
        this.dynamics = new MascotDynamics(stateCount, theta, epochTimes);
        this.maxStep = maxStep;
        this.checkProbabilities = checkProbabilities;

        addModel(treeModel);
        addVariable(tipStates);
        addVariable(theta);
        if (epochTimes != null) {
            addVariable(epochTimes);
        }

        this.eventTapeKnown = false;
        this.coreKnown = false;
        this.likelihoodKnown = false;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            ensureEventTape();
            ensureCore();
            logLikelihood = core.evaluate(eventTape.getEvents(), dynamics.getThetaValues(), false, checkProbabilities).logLikelihood;
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double[] getGradientLogDensity() {
        ensureEventTape();
        ensureCore();
        MascotCore.Result result = core.evaluate(eventTape.getEvents(), dynamics.getThetaValues(), true, checkProbabilities);
        logLikelihood = result.logLikelihood;
        likelihoodKnown = true;
        return result.gradient;
    }

    public Parameter getTheta() {
        return dynamics.getTheta();
    }

    public MascotDynamics getDynamics() {
        return dynamics;
    }

    public int getStateCount() {
        return dynamics.getStateCount();
    }

    public double getMaxStep() {
        return maxStep;
    }

    @Override
    public void makeDirty() {
        eventTapeKnown = false;
        coreKnown = false;
        likelihoodKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            eventTapeKnown = false;
        } else {
            coreKnown = false;
        }
        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == tipStates) {
            eventTapeKnown = false;
        } else if (variable == dynamics.getTheta()) {
            // The event tape is still valid; only the numeric likelihood changes.
        } else if (variable == dynamics.getEpochTimes()) {
            coreKnown = false;
        } else {
            eventTapeKnown = false;
            coreKnown = false;
        }
        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        eventTapeKnown = false;
        coreKnown = false;
    }

    @Override
    protected void acceptState() {
        // Nothing to do.
    }

    @Override
    public String getReport() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    private void ensureEventTape() {
        if (!eventTapeKnown) {
            eventTape = MascotEventTape.fromTree(treeModel, tipStates, dynamics.getStateCount());
            eventTapeKnown = true;
        }
    }

    private void ensureCore() {
        if (!coreKnown) {
            core = new MascotCore(dynamics.getStateCount(), dynamics.getBoundaries(), maxStep);
            coreKnown = true;
        }
    }
}
