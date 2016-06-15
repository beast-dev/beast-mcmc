/*
 * ConditionalCladeProbability.java
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

package dr.evomodel.tree;

import dr.evolution.tree.SimpleTree;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

public class ConditionalCladeProbability extends Likelihood.Abstract {
	
	private ConditionalCladeFrequency ccf;
	private TreeModel treeModel;
	
	public ConditionalCladeProbability(ConditionalCladeFrequency ccf, TreeModel treeModel) {
		super(null);
		this.ccf = ccf;
		this.treeModel = treeModel;
	}

	@Override
	protected double calculateLogLikelihood() {
		
		SimpleTree simTree = new SimpleTree(treeModel);
		
		//System.err.println("tree: " + simTree);
		//System.err.println(ccf.getTreeProbability(simTree));
		
		return ccf.getTreeProbability(simTree);
		
	}
	
	/**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
    }

}
