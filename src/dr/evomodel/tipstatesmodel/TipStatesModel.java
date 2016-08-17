/*
 * TipStatesModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.tipstatesmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public abstract class TipStatesModel extends AbstractModel {

    // an enum which specifies if the model emits tip states or partials
    public enum Type {
        PARTIALS,
        STATES
    };

    /**
     * @param name Model Name
     */
    public TipStatesModel(String name, TaxonList includeTaxa, TaxonList excludeTaxa) {
        super(name);

        this.includeTaxa = includeTaxa;
        this.excludeTaxa = excludeTaxa;
    }

    public final void setTree(Tree tree) {
        this.tree = tree;

        int extNodeCount = tree.getExternalNodeCount();

        excluded = new boolean[extNodeCount];
        if (includeTaxa != null) {
            for (int i = 0; i < extNodeCount; i++) {
                if (includeTaxa.getTaxonIndex(tree.getNodeTaxon(tree.getExternalNode(i))) == -1) {
                    excluded[i] = true;
                }
            }
        }

        if (excludeTaxa != null) {
            for (int i = 0; i < extNodeCount; i++) {
                if (excludeTaxa.getTaxonIndex(tree.getNodeTaxon(tree.getExternalNode(i))) != -1) {
                    excluded[i] = true;
                }
            }

        }

        states = new int[extNodeCount][];

        taxaChanged();
    }

    protected abstract void taxaChanged();

    public final void setStates(PatternList patternList, int sequenceIndex, int nodeIndex, String taxonId) {
        if (this.patternList == null) {
            this.patternList = patternList;
            patternCount = patternList.getPatternCount();
            stateCount = patternList.getDataType().getStateCount();
        } else if (patternList != this.patternList) {
            throw new RuntimeException("The TipStatesModel with id, " + getId() + ", has already been associated with a patternList.");
        }
        if (this.states[nodeIndex] == null) {
            this.states[nodeIndex] = new int[patternCount];
        }

        for (int i = 0; i < patternCount; i++) {
            this.states[nodeIndex][i] = patternList.getPatternState(sequenceIndex, i);
        }

        taxonMap.put(nodeIndex, taxonId);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

    public PatternList getPatternList() {
        return patternList;
    }

    public abstract Type getModelType();

    public abstract void getTipPartials(int nodeIndex, double[] tipPartials);

    public abstract void getTipStates(int nodeIndex, int[] tipStates);

    protected int[][] states;
    protected boolean[] excluded;

    protected int patternCount = 0;
    protected int stateCount;

    protected TaxonList includeTaxa;
    protected TaxonList excludeTaxa;

    protected Tree tree;

    private PatternList patternList = null;

    protected Map<Integer, String> taxonMap = new HashMap<Integer, String>();
}
