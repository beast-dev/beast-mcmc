/*
 * AttributeBranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.math.MathUtils;

import java.util.Iterator;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AttributeBranchRateModel extends AbstractBranchRateModel {

    public AttributeBranchRateModel(final TreeModel treeModel, final String rateAttributeName) {
        super(ATTRIBUTE_BRANCH_RATE_MODEL);

        this.treeModel = treeModel;
        this.rateAttributeName = rateAttributeName;

        addModel(treeModel);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do
    }

    @Override
    protected void storeState() {
        // nothing to do
    }

    @Override
    protected void restoreState() {
        // nothing to do
    }

    @Override
    protected void acceptState() {
        // nothing to do
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        Object value = tree.getNodeAttribute(node, rateAttributeName);
        return Double.parseDouble((String)value);
    }

    @Override
    public String getTraitName() {
        return rateAttributeName;
    }

    public static final String ATTRIBUTE_BRANCH_RATE_MODEL = "attributeBranchRateModel";

    private final TreeModel treeModel;
    private final String rateAttributeName;

}
