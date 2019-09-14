/*
 * PairedParalogBranchModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogBranchModel extends AbstractModel implements BranchModel {

    private final TreeModel tree;
    private final List<SubstitutionModel> substitutionModels;
    private final int[][] assignments;
    private final List<Parameter> timeToDuplicationProportion;
    Map<Integer, Parameter> nodeDuplicationParameterMap = new HashMap<Integer, Parameter>();

    public PairedParalogBranchModel(String name,
                                    List<SubstitutionModel> substitutionModels,
                                    int[][] assignments,
                                    List<Parameter> timeToDuplicationProportions,
                                    TreeModel tree) {
        super(name);
        this.substitutionModels = substitutionModels;
        this.tree = tree;
        this.assignments = assignments;
        this.timeToDuplicationProportion = timeToDuplicationProportions;
        this.nodeDuplicationParameterMap = constructNodeDuplicationTimeMapping();

        for (SubstitutionModel substitutionModel : substitutionModels) {
            addModel(substitutionModel);
        }
        for (Parameter timeProportion : timeToDuplicationProportions) {
            addVariable(timeProportion);
        }
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {
        int[] modelIndices = assignments[branch.getNumber()];
        double[] weights = new double[modelIndices.length];
        double residue = 1.0;
        for (int i = 0; i < weights.length - 1; i++) {
            weights[i] = residue * nodeDuplicationParameterMap.get(branch.getNumber()).getParameterValue(i);
            residue *= 1.0 - weights[i];
        }
        weights[weights.length - 1] = residue;
        return new BranchModel.Mapping() {
            public int[] getOrder() {return modelIndices.clone();}

            @Override
            public double[] getWeights() {
                return weights;
            }
        };
    }

    private Map<Integer, Parameter> constructNodeDuplicationTimeMapping() {
        Map<Integer, Parameter> map = new HashMap<>();
        int j = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            if (assignments[i].length > 1) {
                map.put(i, timeToDuplicationProportion.get(j));
                j++;
            }
        }
        return map;
    }

    public static int[][] parseAssignmentString(String[] assignmentString) {
        int[][] assignment = new int[assignmentString.length][];
        List<Integer> modelList = new ArrayList<>();
        for (int i = 0; i < assignmentString.length; i++) {
            String thisAssignment = assignmentString[i];
            int start = 0; int end = 0;
            modelList.clear();
            while(end < thisAssignment.length()) {
                Character character = thisAssignment.charAt(end);
                if( character.equals('+')) {
                    final int modelIndex = Integer.valueOf(thisAssignment.substring(start, end));
                    modelList.add(modelIndex);
                    start = end + 1;
                }
                end++;
            }
            final int modelIndex = Integer.valueOf(thisAssignment.substring(start, end));
            modelList.add(modelIndex);
            assignment[i] = new int[modelList.size()];
            for (int j = 0; j < modelList.size(); j++) {
                assignment[i][j] = modelList.get(j);
            }
        }
        return assignment;
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModels;
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return substitutionModels.get(assignments[tree.getRoot().getNumber()][0]);
    }

    @Override
    public FrequencyModel getRootFrequencyModel() {
        return getRootSubstitutionModel().getFrequencyModel();
    }

    @Override
    public boolean requiresMatrixConvolution() {
        return true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
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
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }
}
