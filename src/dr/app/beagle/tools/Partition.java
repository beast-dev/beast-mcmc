/*
 * Partition.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
package dr.app.beagle.tools;

import java.util.HashMap;
import java.util.Map;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.bss.Utils;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class Partition {

	public int from;
	public int to;
	public int every;

	private Integer partitionNumber; 
	private BranchModel branchModel;
	private TreeModel treeModel;
	private GammaSiteRateModel siteModel;
	private BranchRateModel branchRateModel;
	private FrequencyModel freqModel;
	
	public boolean hasAncestralSequence = false;
	public Sequence ancestralSequence = null;
	
	private Map<Taxon,int[]> sequenceList; 
	
	public Partition(TreeModel treeModel, //
			BranchModel branchModel, //
			GammaSiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int from, //
			int to, //
			int every //
	) {

		this.setTreeModel(treeModel);
		this.setSiteModel(siteModel);
		this.setFreqModel(freqModel);
		this.branchModel = branchModel;
		this.setBranchRateModel(branchRateModel);

		this.from = from;
		this.to = to;
		this.every = every;
		
		sequenceList = new HashMap<Taxon,int[]>();
		
	}// END: Constructor

	public void setAncestralSequence(Sequence ancestralSequence) {
		this.ancestralSequence = ancestralSequence;
		this.hasAncestralSequence = true;
	}// END: setAncestralSequence

	public int getPartitionSiteCount() {
		return ((to - from) / every) + 1;
	}// END: getPartitionSiteCount

	public BranchModel getBranchModel() {
		return this.branchModel;
	}// END: getBranchModelic 
	
	public void printSequences() {
		
		System.out.println("partition "+ partitionNumber);
		Utils.printMap(getSequenceList());
	}

	public void setPartitionNumber(Integer partitionNumber) {
		this.partitionNumber = partitionNumber;
	}

	public Integer getPartitionNumber() {
		return partitionNumber;
	}

	public TreeModel getTreeModel() {
		return treeModel;
	}

	public void setTreeModel(TreeModel treeModel) {
		this.treeModel = treeModel;
	}

	public GammaSiteRateModel getSiteModel() {
		return siteModel;
	}

	public void setSiteModel(GammaSiteRateModel siteModel) {
		this.siteModel = siteModel;
	}

	public FrequencyModel getFreqModel() {
		return freqModel;
	}

	public void setFreqModel(FrequencyModel freqModel) {
		this.freqModel = freqModel;
	}

	public Sequence getAncestralSequence() {
		return ancestralSequence;
	}

	public Map<Taxon,int[]> getSequenceList() {
		return sequenceList;
	}

	public BranchRateModel getBranchRateModel() {
		return branchRateModel;
	}

	public void setBranchRateModel(BranchRateModel branchRateModel) {
		this.branchRateModel = branchRateModel;
	}
	
}// END: class
