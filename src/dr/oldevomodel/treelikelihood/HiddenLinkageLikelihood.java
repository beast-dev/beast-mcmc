/*
 * HiddenLinkageLikelihood.java
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.LinkedGroup;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

import java.util.ArrayList;
import java.util.Set;

/**
 * @author Aaron Darling
 */
@Deprecated // Switching to BEAGLE
public class HiddenLinkageLikelihood extends AbstractModelLikelihood
{

	TreeModel tree;
	HiddenLinkageModel hlm;
	

	public HiddenLinkageLikelihood(HiddenLinkageModel hlm, TreeModel tree)
	{
		super("HiddenLinkageLikelihood");
		this.hlm = hlm;
		this.tree = tree;
	}



	protected void acceptState() {
		// nothing to do
	}



	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// nothing to do
	}



	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		// nothing to do
	}



	protected void restoreState() {
		// nothing to do
	}



	protected void storeState() {
		// nothing to do
	}



	public double getLogLikelihood() {
		double logL = 0;
		// first check whether the reads are linked together
		// according to the linkage constraints provided as data in the
		// input file
		ArrayList<LinkedGroup> linkedGroups = hlm.getData().getConstraints();
		for( LinkedGroup lg : linkedGroups ){
			TaxonList tl = lg.getLinkedReads();
			int found = 0;
			for( int l=0; l<hlm.getLinkageGroupCount(); l++ ){
				Set<Taxon> group = hlm.getGroup(l);
				for(int i=0; i<tl.getTaxonCount(); i++){
					if(group.contains(tl.getTaxon(i)))
						found++;
				}
				if(found==tl.getTaxonCount()){
					logL += Math.log(lg.getLinkageProbability());
					break;
				}else if(found>0){
					logL += Math.log(1.0-lg.getLinkageProbability());
					break;
				}
			}
		}

		// then check whether the topology of reference taxa is consistent with the fixed 
		// reference tree topology
		if(hlm.getData().getFixedReferenceTree()){
			logL += Double.NEGATIVE_INFINITY;	// not yet implemented.
		}
		return logL;
	}



	public Model getModel() {
		return this;
	}


	
	public void makeDirty() {
		// nothing to do
	}

}
