/*
 * BranchRateGradientForDiscreteTrait.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.*;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BranchRateGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient implements Citable {

    public BranchRateGradientForDiscreteTrait(String traitName,
                                              TreeDataLikelihood treeDataLikelihood,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              Parameter rateParameter,
                                              boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
    }

    @Override
    protected double getChainGradient(Tree tree, NodeRef node) {
        return branchRateModel.getBranchRateDifferential(tree, node) * tree.getBranchLength(node);
    }

    @Override
    protected double getChainSecondDerivative(Tree tree, NodeRef node) {
        return branchRateModel.getBranchRateSecondDifferential(tree, node) * tree.getBranchLength(node);
    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using linear-time branch-rate differential calculations";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author( "Z", "Zhang"),
                    new Author("A", "Holbrook"),
                    new Author("A", "Nishimura"),
                    new Author("G", "Beale"),
                    new Author("A", "Rambaut"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Gradients do grow on trees: a linear-time O(N)-dimensional gradient for statistical phylogenetics",
            "Annals of Applied Statistics",
            Citation.Status.IN_SUBMISSION);
}
