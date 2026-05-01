package dr.evomodel.treedatalikelihood.continuous.canonical.contribution;

import dr.evomodel.treedatalikelihood.continuous.canonical.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.traversal.*;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.*;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;

interface BranchContributionStrategy {

    boolean contributes();

    boolean fillContribution(int childIndex,
                             CanonicalGaussianTransition transition,
                             CanonicalBranchContributionAssembler assembler,
                             BranchGradientWorkspace workspace);
}
