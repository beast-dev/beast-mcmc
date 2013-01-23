/*
 * CTMCScalePrior.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.MSSD;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.GammaFunction;

/**
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * @author Marc A. Suchard
 * 
 * Date: Aug 22, 2008
 * Time: 3:26:57 PM
 */
public class CTMCScalePrior extends AbstractModelLikelihood {
    final private Parameter ctmcScale;
    final private TreeModel treeModel;
    private double treeLength;
    private boolean treeLengthKnown;

    private static final double logGammaOneHalf = GammaFunction.lnGamma(0.5);

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel) {
        super(name);
        this.ctmcScale = ctmcScale;
        this.treeModel = treeModel;
        addModel(treeModel);
        treeLengthKnown = false;
    }

    private void updateTreeLength() {
        treeLength = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            treeLengthKnown = false;
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
        treeLengthKnown = false;
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        double ab = ctmcScale.getParameterValue(0);
//        if (!treeLengthKnown) {
//            updateTreeLength();
//            treeLengthKnown = true;
//        }
//        double totalTreeTime = treeLength;
        double totalTreeTime = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
        double logNormalization = 0.5 * Math.log(totalTreeTime) - logGammaOneHalf;
        return logNormalization - 0.5 * Math.log(ab) - ab * totalTreeTime; // TODO Change to treeLength and confirm results
    }

    public void makeDirty() {
        treeLengthKnown = false;
    }
}
