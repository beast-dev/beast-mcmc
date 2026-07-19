/*
 * AbstractStructuredCoalescentLikelihood.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;

/**
 * Scaffolding shared by BASTA's {@code BastaLikelihood} and MASCOT's {@code
 * MascotLikelihood}: both own a {@link BestSignalsFromBigFastTreeIntervals}
 * as the single source of truth for "this tree's coalescent/sample events in
 * time order," and an optional {@link BranchRateModel} clock. Everything
 * downstream of that -- interval traversal, transition-matrix/ODE evaluation,
 * gradients -- is genuinely different between the two engines (matrix
 * exponentiation vs. RK4 integration) and stays in each subclass.
 */
public abstract class AbstractStructuredCoalescentLikelihood extends AbstractModelLikelihood {

    protected final BestSignalsFromBigFastTreeIntervals treeIntervals;
    // Nullable: not every structured-coalescent likelihood requires a clock
    // (MASCOT allows null, meaning no rate scaling of the migration process).
    protected final BranchRateModel branchRateModel;

    protected AbstractStructuredCoalescentLikelihood(String name,
                                                       BestSignalsFromBigFastTreeIntervals treeIntervals,
                                                       BranchRateModel branchRateModel) {
        super(name);
        this.treeIntervals = treeIntervals;
        addModel(treeIntervals);
        this.branchRateModel = branchRateModel;
        if (branchRateModel != null) {
            addModel(branchRateModel);
        }
    }

    @Override
    public final Model getModel() {
        return this;
    }
}
