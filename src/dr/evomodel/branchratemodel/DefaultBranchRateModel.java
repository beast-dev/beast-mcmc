/*
 * DefaultBranchRateModel.java
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

import dr.evolution.tree.*;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Variable;

/**
 * @author Andrew Rambaut
 * @version $Id: DefaultBranchRateModel.java,v 1.4 2005/05/24 20:25:57 rambaut Exp $
 */
public final class DefaultBranchRateModel implements BranchRateModel {
    public double getBranchRate(Tree tree, NodeRef node) {
        return 1.0;
    }

    public void addModelListener(ModelListener listener) {
        // nothing to do
    }

    public void removeModelListener(ModelListener listener) {
        // nothing to do
    }

    public void storeModelState() {
        // nothing to do
    }

    public void restoreModelState() {
        // nothing to do
    }

    public void acceptModelState() {
        // nothing to do
    }

    public boolean isValidState() {
        return true;
    }

    public int getModelCount() {
        return 0;
    }

    public Model getModel(int i) {
        return null;
    }

    public int getVariableCount() {
        return 0;
    }

    public Variable getVariable(int i) {
        return null;
    }

    public String getModelName() {
        return null;
    }

    public String getId() {
        return null;
    }

    public void setId(String id) {
        // nothing to do
    }

    public boolean isUsed() {
        return false;
    }

    public String getTraitName() {
        return RATE;
    }

    public Intent getIntent() {
        return Intent.BRANCH;
    }

    public boolean getLoggable() {
        return true;
    }

    public TreeTrait getTreeTrait(final String key) {
        return this;
    }

    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[] { this };
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public int getDimension() {
        return 1;
    }

    public Double getTrait(final Tree tree, final NodeRef node) {
        return getBranchRate(tree, node);
    }

    public String getTraitString(final Tree tree, final NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }
}
