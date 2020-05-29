/*
 * RzhetskyNeiBranchLengthsTree.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import java.util.*;

/**
 * A tree in which the branch lengths are calculated from a distance matrix
 * using the methods of Rzhetsky & Nei (1993) MBE.
 *
 * @author Andrew Rambaut
 */
public class RzhetskyNeiBranchLengthsTree extends SimpleTree {

	/**
	 * constructor
	 */
	public RzhetskyNeiBranchLengthsTree(Tree sourceTree, DistanceMatrix distanceMatrix) {

		super(sourceTree);

		this.distanceMatrix = distanceMatrix;

		Map<NodeRef, Set<Integer>> taxonSetMap = new HashMap<>();

		allTaxonSet = new HashSet<>(getTaxonSets(this, getRoot(), taxonSetMap));
	}

	private Set<Integer> getTaxonSets(Tree tree, NodeRef node, Map<NodeRef, Set<Integer>> taxonSetMap) {

		Set<Integer> taxonSet = new HashSet<>();
		if (tree.isExternal(node)) {
			taxonSet.add(node.getNumber());
		} else {
			assert tree.getChildCount(node) == 2 : "Must be a strictly bifurcating tree";

			for (int i = 0; i < tree.getChildCount(node); i++) {
				taxonSet.addAll(getTaxonSets(tree, getChild(node, i), taxonSetMap));
			}
		}

		taxonSetMap.put(node, taxonSet);

		return taxonSet;
	}

	private void calculateBranchLengths(Tree tree, NodeRef node, NodeRef sibling, Map<NodeRef, Set<Integer>> taxonSetMap) {

		double length;

		if (tree.isExternal(node)) {
			Set<Integer> taxonSetC = taxonSetMap.get(node);
			Set<Integer> taxonSetB = taxonSetMap.get(sibling);

			Set<Integer> taxonSetA = new HashSet<>(allTaxonSet);
			taxonSetA.removeAll(taxonSetC);
			taxonSetA.removeAll(taxonSetB);

			double nA = taxonSetA.size();
			double nB = taxonSetB.size();

			double dCA = getSumOfDistances(taxonSetC, taxonSetA);
			double dCB = getSumOfDistances(taxonSetC, taxonSetB);
			double dAB = getSumOfDistances(taxonSetA, taxonSetB);

			// Equation 4 from R&N1993
			length = 0.5 * ((dCA / nA) + (dCB / nB) - (dAB / (nA * nB)));

		} else {

			NodeRef child1 = getChild(node, 0);
			NodeRef child2 = getChild(node, 1);

			calculateBranchLengths(tree, child1, child2, taxonSetMap);

			calculateBranchLengths(tree, child2, child1, taxonSetMap);

			Set<Integer> taxonSetC = taxonSetMap.get(child1);
			Set<Integer> taxonSetD = taxonSetMap.get(child2);

			Set<Integer> taxonSetB = taxonSetMap.get(sibling);

			Set<Integer> taxonSetA = new HashSet<>(allTaxonSet);
			taxonSetA.removeAll(taxonSetC);
			taxonSetA.removeAll(taxonSetD);
			taxonSetA.removeAll(taxonSetB);

			double nA = taxonSetA.size();
			double nB = taxonSetB.size();
			double nC = taxonSetC.size();
			double nD = taxonSetD.size();

			// Equation 3 from R&N1993
			double gamma = (nB * nC + nA * nD) / ((nA + nB)*(nC + nD));

			double dAC = getSumOfDistances(taxonSetA, taxonSetC);
			double dBD = getSumOfDistances(taxonSetB, taxonSetD);
			double dBC = getSumOfDistances(taxonSetB, taxonSetC);
			double dAD = getSumOfDistances(taxonSetA, taxonSetD);
			double dAB = getSumOfDistances(taxonSetA, taxonSetB);
			double dCD = getSumOfDistances(taxonSetC, taxonSetD);

			// Equation 2 from R&N1993
			length = 0.5 * (
					gamma * ((dAC / nA * nC) + (dBD / nB * nD)) +
							(1.0 - gamma) * ((dBC / nB * nC) + (dAD / nA * nD)) -
							(dAB / nA * nB) -
							(dCD / nC * nD)
					);
		}

		setBranchLength(node, length);
	}

	private double getSumOfDistances(Set<Integer> taxonSet1, Set<Integer> taxonSet2) {
		double sum = 0.0;

		for (int i : taxonSet1) {
			for (int j : taxonSet2) {
				sum += distanceMatrix.getElement(i, j);
			}
		}

		return sum;
	}

	private final DistanceMatrix distanceMatrix;
	private final Set<Integer> allTaxonSet;
}
