/*
 * BranchModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.inference.model.Model;

import java.util.List;

/**
 * This is an interface which provides a mapping of substitution models to branches.
 *
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
public interface BranchModel extends Model  {
    /**
     * Returns a mapping of substitution models to the given branch. The Mapping
     * contains a list of substitution models in order from tipward to rootward
     * and a set of relative weights for each (may be times or proportions).
     *
     * @param branch the branch
     * @return a Mapping object
     */
    Mapping getBranchModelMapping(final NodeRef branch);

    /**
     * Gets the list of substitution models in order they will be referred to
     * by the indices returned by the mappings.
     * @return the list of substitution models
     */
    List<SubstitutionModel> getSubstitutionModels();

    /**
     * Gets the substitution model that will be applied at the root.
     * @return the substitution model
     */
    SubstitutionModel getRootSubstitutionModel();

    /**
     * Gets the frequency model that will be applied at the root.
     * @return the substitution model
     */
    FrequencyModel getRootFrequencyModel();

    /**
     * Is this model going to require convolution of matrices along any branches (essentially
     * are the mappings ever going to return more than one substitution model.
     * @return does it?
     */
    boolean requiresMatrixConvolution();

    interface Mapping {
        int[] getOrder();
        double[] getWeights();
    }

    Mapping DEFAULT = new Mapping() {
        public int[] getOrder() {
            return new int[] { 0 };
        }

        public double[] getWeights() {
            return new double[] { 1.0 };
        }
    };
}
