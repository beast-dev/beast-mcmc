/*
 * ContractedRatesProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.substmodel;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Bounds;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * One rate per internal node of a fixed multifurcating character tree, given in the post-order of the Newick as
 * written (children left to right). The model runs on the binary skeleton: a kept skeleton node reads the rate of its
 * clade, a collapsed one has rate 0. The rates are normalised to one expected substitution per unit time,
 * mu = sum_m r_m (pi_m - sum_{i in m} pi_i^2 / pi_m), and bounded below by the valid-generator cone of the
 * multifurcating tree (docs/generator-constraints), so negative rates are allowed where the generator stays valid.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ContractedRatesProvider extends AbstractModel implements HierarchicalSubstitutionModel.HierarchicalRateProvider {

    private final TreeModel skeleton;
    private final Parameter rates;
    private final FrequencyModel frequency;
    private final CharacterTreeContraction contraction;
    private final int tips;
    private final int[] slotOfInternal;
    private final int[][] statesOfSlot;
    private final int[] parentSlot;
    private final int[][] childSlots;
    private final double[] lastRates;
    private final double[] lastFrequencies;
    private double normalization = Double.NaN;

    public ContractedRatesProvider(Parameter rates, CharacterTreeContraction contraction, FrequencyModel frequency,
                                   int[] slotOfInternal) {
        super("contractedRates");
        if (contraction.getKeep() != null) {
            throw new IllegalArgumentException("one rate per internal node of a multifurcating tree needs a fixed tree");
        }
        this.skeleton = contraction.getTree();
        this.rates = rates;
        this.frequency = frequency;
        this.contraction = contraction;
        this.tips = skeleton.getExternalNodeCount();
        if (slotOfInternal.length != skeleton.getInternalNodeCount()) {
            throw new IllegalArgumentException("one slot per internal node of the skeleton is needed");
        }
        this.slotOfInternal = slotOfInternal.clone();
        final int slots = contraction.getContractedInternalNodeCount();
        if (rates.getDimension() != slots) {
            throw new IllegalArgumentException("ratesParameter has dimension " + rates.getDimension()
                    + " but the character tree has " + slots + " internal nodes");
        }
        this.statesOfSlot = new int[slots][];
        this.parentSlot = new int[slots];
        this.childSlots = new int[slots][];
        for (int i = 0; i < slotOfInternal.length; i++) {
            final int node = tips + i;
            if (contraction.isCollapsed(node) != (slotOfInternal[i] < 0)) {
                throw new IllegalArgumentException("the slots must be the kept nodes of the character tree");
            }
            final int slot = slotOfInternal[i];
            if (slot < 0) {
                continue;
            }
            if (slot >= slots || statesOfSlot[slot] != null) {
                throw new IllegalArgumentException("slot " + slot + " is not a distinct internal node of the character tree");
            }
            statesOfSlot[slot] = contraction.sortedStatesBelow(node);
            final int parent = contraction.contractedParent(node);
            parentSlot[slot] = parent < 0 ? -1 : slotOfInternal[parent - tips];
            final List<Integer> children = new ArrayList<Integer>();
            for (int child : contraction.contractedChildren(node)) {
                if (child >= tips) {
                    children.add(slotOfInternal[child - tips]);
                }
            }
            childSlots[slot] = new int[children.size()];
            for (int c = 0; c < children.size(); c++) {
                childSlots[slot][c] = children.get(c);
            }
        }
        this.lastRates = new double[slots];
        this.lastFrequencies = new double[frequency.getDataType().getStateCount()];
        rates.addBounds(new GeneratorBounds());
        addVariable(rates);
        addModel(frequency);
        addModel(skeleton);
    }

    /**
     * The rate of each internal node of the multifurcating tree, in the post-order of its Newick.
     */
    public Parameter getRatesParameter() {
        return rates;
    }

    public int getSlotOfNode(int nodeNumber) {
        return nodeNumber < tips ? -1 : slotOfInternal[nodeNumber - tips];
    }

    public double getNormalization() {
        boolean changed = Double.isNaN(normalization);
        for (int i = 0; i < lastRates.length; i++) {
            final double r = rates.getParameterValue(i);
            if (r != lastRates[i]) {
                lastRates[i] = r;
                changed = true;
            }
        }
        for (int i = 0; i < lastFrequencies.length; i++) {
            final double f = frequency.getFrequency(i);
            if (f != lastFrequencies[i]) {
                lastFrequencies[i] = f;
                changed = true;
            }
        }
        if (changed) {
            double mu = 0.0;
            for (int node : contraction.getKeptInternalNodes()) {
                final int slot = slotOfInternal[node - tips];
                double piClade = 0.0, sumSquares = 0.0;
                for (int state : statesOfSlot[slot]) {
                    piClade += lastFrequencies[state];
                    sumSquares += lastFrequencies[state] * lastFrequencies[state];
                }
                mu += lastRates[slot] * (piClade - sumSquares / piClade);
            }
            normalization = mu;
        }
        return normalization;
    }

    @Override
    public double getRate(TreeModel tree, int node) {
        if (contraction.isCollapsed(node)) {
            return 0.0;
        }
        return rates.getParameterValue(slotOfInternal[node - tips]) / getNormalization();
    }

    @Override
    public TreeModel getTreeModel() {
        return skeleton;
    }

    @Override
    public FrequencyModel getFrequencyModel() {
        return frequency;
    }

    @Override
    public CharacterTreeContraction getContraction() {
        return contraction;
    }

    @Override
    public void addToModelListener(AbstractModel listener) {
        listener.addModel(this);
    }

    /**
     * The exact lower limits of the rates on the multifurcating tree, as GeneratorBounds of the Default provider:
     * with the other rates fixed, r_m must keep sigma(n) >= 0 at m and at every internal node below it.
     */
    private final class GeneratorBounds implements Bounds<Double> {

        private double cladeFrequency(int slot) {
            double sum = 0.0;
            for (int state : statesOfSlot[slot]) {
                sum += frequency.getFrequency(state);
            }
            return sum;
        }

        private double term(int slot) {
            return rates.getParameterValue(slot) / cladeFrequency(slot);
        }

        private double minimumBelow(int slot, double above) {
            double minimum = Double.POSITIVE_INFINITY;
            for (int c : childSlots[slot]) {
                final double sum = above + term(c);
                minimum = Math.min(minimum, Math.min(sum, minimumBelow(c, sum)));
            }
            return minimum;
        }

        @Override
        public Double getUpperLimit(int dimension) {
            return Double.POSITIVE_INFINITY;
        }

        @Override
        public Double getLowerLimit(int dimension) {
            double sigmaParent = 0.0;
            for (int k = parentSlot[dimension]; k >= 0; k = parentSlot[k]) {
                sigmaParent += term(k);
            }
            final double below = Math.min(0.0, minimumBelow(dimension, 0.0));
            return -cladeFrequency(dimension) * (sigmaParent + below);
        }

        @Override
        public int getBoundsDimension() {
            return rates.getDimension();
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == skeleton) {
            if (HierarchicalSubstitutionModel.changesTopology(object)) {
                throw new IllegalStateException("the topology of the character tree " + skeleton.getId()
                        + " cannot change with one rate per internal node of a fixed multifurcating tree");
            }
            return;
        }
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public String toString() {
        return "ContractedRatesProvider" + Arrays.toString(rates.getParameterValues());
    }
}
