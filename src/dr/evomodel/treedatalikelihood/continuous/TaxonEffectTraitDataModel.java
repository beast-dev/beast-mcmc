/*
 * IntegratedProcessTraitDataModel.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class TaxonEffectTraitDataModel extends
        ContinuousTraitDataModel implements ContinuousTraitPartialsProvider {

    private final EffectMap map;
    private final Parameter effects;
    private final Tree tree;

    public TaxonEffectTraitDataModel(String name,
                                     Tree tree,
                                     CompoundParameter parameter,
                                     EffectMap map,
                                     boolean[] missingIndicators,
                                     boolean useMissingIndices,
                                     final int dimTrait,
                                     PrecisionType precisionType) {
        super(name, parameter, missingIndicators, useMissingIndices, dimTrait, precisionType);

        this.tree = tree;
        this.map = map;
        this.effects = map.parameter;

        if (parameter.getDimension() != effects.getDimension()) {
            throw new IllegalArgumentException("Unequal effects dimension");
        }

        if (dimTrait != 1) {
            throw new IllegalArgumentException("Not yet implemented");
        }

        addVariable(effects);
    }

    public EffectMap getMap() {
        return map;
    }

    public Tree getTree() {
        return tree;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);

        int index = map.getEffectIndex(taxonIndex);

//        if (index == 1) {
//            System.err.println("one");
//        }
        double effect = effects.getParameterValue(index);
        partial[0] -= effect;

        return partial;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.handleVariableChangedEvent(variable, index, type);
        if (variable == effects) {
            if (type == Parameter.ChangeType.VALUE_CHANGED) {
                fireModelChanged(this, getTaxonIndex(index));
            } else if (type == Parameter.ChangeType.ALL_VALUES_CHANGED) {
                fireModelChanged(this);
            } else {
                throw new RuntimeException("Unhandled parameter change type");
            }
        }
    }

    @Override
    public void updateTipDataGradient(DenseMatrix64F precision, DenseMatrix64F variance, NodeRef node,
                                      int offset, int dimGradient) {
        throw new RuntimeException("not yet implemented");
    }

    @Override
    public boolean needToUpdateTipDataGradient(int offset, int dimGradient) {
        throw new RuntimeException("not yet implemented");
    }

    public static class EffectMap {

        private final Parameter parameter;
        private final Map<String,Integer> taxonNameToIndex;
        private final int[] map;
        private final int dim;

        public EffectMap(Tree tree, Parameter parameter, int dim) {

            if (tree.getExternalNodeCount() != parameter.getDimension() * dim) {
                throw new IllegalArgumentException("Invalid effect dimension");
            }

            this.parameter = parameter;
            this.dim = dim;
            this.taxonNameToIndex = new HashMap<>();

            String base = parameter.getParameterName();

            String[] names = new String[tree.getExternalNodeCount()];
            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                String taxonName = tree.getTaxonId(i);
                taxonNameToIndex.put(taxonName, i);
                for (int j = 0; j < dim; ++j) {
                    names[dim * i + j] = makeName(base, taxonName, j + 1);
                }
            }

            parameter.setDimensionNames(names);

            this.map = makeMap(tree, taxonNameToIndex);
        }

        public EffectMap(Tree tree, EffectMap original) {

            this.parameter = original.parameter;
            this.taxonNameToIndex = original.taxonNameToIndex;
            this.dim = original.dim;

            this.map = makeMap(tree, taxonNameToIndex);
        }

        public int getEffectIndex(int taxonId) {
            return map[taxonId];
        }

        public Parameter getEffects() { return parameter; }

        private int[] makeMap(Tree tree, Map<String,Integer> taxonNameToIndex) {
            int[] map = new int[tree.getExternalNodeCount()];
            Arrays.fill(map, -1);

            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                String taxonId = tree.getTaxonId(i);
                Integer value = taxonNameToIndex.get(taxonId);
                if (value == null) {
                    throw new IllegalArgumentException("Unable to find taxon '" + taxonId + "'");
                }
                map[i] = value;
            }

            return map;
        }

        private static String makeName(String base, String taxonName, int dimension) {
            return base + "." + taxonName + "." + dimension;
        }
    }
}