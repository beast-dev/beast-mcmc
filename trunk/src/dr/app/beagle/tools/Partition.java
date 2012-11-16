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

package dr.app.beagle.tools;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.evolution.sequence.Sequence;
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

	private BranchModel branchModel;
	
	public TreeModel treeModel;
	public GammaSiteRateModel siteModel;
	public BranchRateModel branchRateModel;
	public FrequencyModel freqModel;

	public boolean hasAncestralSequence = false;
	public Sequence ancestralSequence = null;

	// public static final ThreadLocal threadLocal = new ThreadLocal();

	public Partition(TreeModel treeModel, //
			BranchModel branchModel, GammaSiteRateModel siteModel, //
			BranchRateModel branchRateModel, //
			FrequencyModel freqModel, //
			int from, //
			int to, //
			int every //
	) {

		this.treeModel = treeModel;
		this.siteModel = siteModel;
		this.freqModel = freqModel;
		this.branchModel = branchModel;
		this.branchRateModel = branchRateModel;

		this.from = from;
		this.to = to;
		this.every = every;

	}// END: Constructor

	public void setAncestralSequence(Sequence ancestralSequence) {
		this.ancestralSequence = ancestralSequence;
		this.hasAncestralSequence = true;
	}// END: setAncestralSequence

	public int getPartitionSiteCount() {
		int partitionSiteCount = ((to - from) / every) + 1;
		return partitionSiteCount;
	}// END: getPartitionSiteCount

	public BranchModel getBranchModel() {
		return this.branchModel;
	}// END: getBranchModel

}// END: class
