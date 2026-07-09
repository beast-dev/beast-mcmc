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
    private boolean gradientKnown;

    private double logLikelihood;
    private double[] gradient;

    // Snapshots taken by storeState() and restored by restoreState() on a
    // rejected proposal. Restoring the actual eventTape/core references (rather
    // than just invalidating eventTapeKnown/coreKnown, as an earlier version of
    // this class did) means a rejected proposal that never touched the tree or
    // epoch times does not force a full event-tape/core rebuild on the next
    // evaluation. Both MascotEventTape and MascotCore are effectively immutable
    // from this wrapper's point of view once built (MascotCore's internal
    // per-epoch rate cache is fully recomputed from theta at the start of every
    // evaluate() call regardless), so sharing the stored reference back in is safe.
    private MascotEventTape storedEventTape;
    private MascotCore storedCore;
    private boolean storedEventTapeKnown;
    private boolean storedCoreKnown;
    private boolean storedLikelihoodKnown;
    private boolean storedGradientKnown;
    private double storedLogLikelihood;
    private double[] storedGradient;

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
        this.gradientKnown = false;
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
            logLikelihood = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(), false, checkProbabilities, false).logLikelihood;
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double[] getGradientLogDensity() {
        if (!gradientKnown) {
            ensureEventTape();
            ensureCore();
            MascotCore.Result result = core.evaluate(eventTape.getPreparedEvents(), dynamics.getThetaValues(), true, checkProbabilities, false);
            logLikelihood = result.logLikelihood;
            gradient = result.gradient;
            likelihoodKnown = true;
            gradientKnown = true;
        }
        // Clone on every return, not just on a cache hit: the stored gradient array
        // is reused across calls, so a caller that mutated the returned array in
        // place would otherwise silently corrupt the cache for the next caller.
        return gradient.clone();
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
        gradientKnown = false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            eventTapeKnown = false;
        } else {
            coreKnown = false;
        }
        likelihoodKnown = false;
        gradientKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == tipStates) {
            eventTapeKnown = false;
        } else if (variable == dynamics.getTheta()) {
            // The event tape is still valid; only the numeric likelihood and gradient change.
        } else if (variable == dynamics.getEpochTimes()) {
            coreKnown = false;
        } else {
            eventTapeKnown = false;
            coreKnown = false;
        }
        likelihoodKnown = false;
        gradientKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedGradient = gradient;
        storedEventTape = eventTape;
        storedCore = core;
        storedEventTapeKnown = eventTapeKnown;
        storedCoreKnown = coreKnown;
        storedLikelihoodKnown = likelihoodKnown;
        storedGradientKnown = gradientKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        gradient = storedGradient;
        eventTape = storedEventTape;
        core = storedCore;
        eventTapeKnown = storedEventTapeKnown;
        coreKnown = storedCoreKnown;
        likelihoodKnown = storedLikelihoodKnown;
        gradientKnown = storedGradientKnown;
    }

    @Override
    protected void acceptState() {
        // Nothing to do: current fields already reflect the accepted state.
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
