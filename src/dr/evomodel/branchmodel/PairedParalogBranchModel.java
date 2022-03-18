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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.DuplicationTreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogBranchModel extends AbstractModel implements BranchModel {

    private final DuplicationTreeModel tree;
    private final SubstitutionModel baseModel;
    private final SubstitutionModel geneConversionModel;
    private final Taxa postDuplicationTaxa;
    private final List<SubstitutionModel> substitutionModels;
    private final Map<Integer, SubstitutionModel> branchModelMap;

    public PairedParalogBranchModel(String name,
                                    SubstitutionModel baseModel,
                                    SubstitutionModel geneConversionModel,
                                    DuplicationTreeModel tree,
                                    Taxa postDuplicationTaxa) {
        super(name);

        this.tree = tree;
        this.baseModel = baseModel;
        this.geneConversionModel = geneConversionModel;
        this.postDuplicationTaxa = postDuplicationTaxa;
        this.substitutionModels = new ArrayList<>();
        substitutionModels.add(baseModel);
        substitutionModels.add(geneConversionModel);
        this.branchModelMap = new HashMap<>();

        buildBranchModelMapping();

        addModel(baseModel);
        addModel(geneConversionModel);

    }

    private void buildBranchModelMapping() {
        Set<String> leafNodes = new HashSet<>();

        for (int i = 0; i < postDuplicationTaxa.getTaxonCount(); i++) {
            leafNodes.add(postDuplicationTaxa.getTaxon(i).getId());
        }

        NodeRef duplicationNode = TreeUtils.getCommonAncestorNode(tree, leafNodes);

        buildBranchModelMappingRecursively(tree, tree.getRoot(), duplicationNode, branchModelMap, duplicationNode == tree.getRoot());
    }

    private void buildBranchModelMappingRecursively(Tree tree, NodeRef node, NodeRef duplicationNode,
                                                    Map<Integer, SubstitutionModel> map, boolean postDuplication) {
        if (postDuplication) {
            map.put(node.getNumber(), geneConversionModel);
        } else {
            map.put(node.getNumber(), baseModel);
        }
        if (node == duplicationNode) {
            postDuplication = true;
        }
        for (int i = 0; i < tree.getChildCount(node); i++) {
            buildBranchModelMappingRecursively(tree, tree.getChild(node, i), duplicationNode, map, postDuplication);
        }
    }

    private int getSubstitutionModelIndex(SubstitutionModel substitutionModel) {
        return substitutionModel == baseModel ? 0 : 1;
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {
        return new BranchModel.Mapping() {
            public int[] getOrder() {return new int[]{getSubstitutionModelIndex(branchModelMap.get(branch.getNumber()))};}

            @Override
            public double[] getWeights() {
                return new double[]{1.0};
            }
        };
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModels;
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return branchModelMap.get(tree.getRoot().getNumber());
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
        fireModelChanged(object, index);
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
