/*
 * UPGMATree.java
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

package dr.evolution.tree;

import dr.evolution.distance.DistanceMatrix;

/**
 * constructs a UPGMA tree from pairwise distances
 *
 * @version $Id: UPGMATree.java,v 1.14 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class UPGMATree extends ClusteringTree {

	/**
	 * constructor UPGMA tree
	 *
	 * @param distanceMatrix distance matrix
	 */
	public UPGMATree(DistanceMatrix distanceMatrix) {
	
		super(distanceMatrix, 2);
	}

	//
	// Protected and Private stuff
	//

	protected void findNextPair() {
	
		besti = 0;
		bestj = 1;
		double dmin = getDist(0, 1);
		for (int i = 0; i < numClusters-1; i++) {
		
			for (int j = i+1; j < numClusters; j++) {
			
				if (getDist(i, j) < dmin) {
				
					dmin = getDist(i, j);
					besti = i;
					bestj = j;
				}
			}
		}
		abi = alias[besti];
		abj = alias[bestj];
	}

	protected double newNodeHeight() {
		return getDist(besti, bestj) / 2.0;
	}

	/**
	 * compute updated distance between the new cluster (i,j)
	 * to any other cluster k
	 */
	protected double updatedDistance(int i, int j, int k)
	{
		int ai = alias[i];
		int aj = alias[j];
		
		double tipSum = (double) (tipCount[ai] + tipCount[aj]);
		
		return 	(((double)tipCount[ai]) / tipSum) * getDist(k, i) +
				(((double)tipCount[aj]) / tipSum) * getDist(k, j);
	}
}
