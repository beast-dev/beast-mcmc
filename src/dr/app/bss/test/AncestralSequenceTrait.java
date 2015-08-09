/*
 * AncestralSequenceTrait.java
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

package dr.app.bss.test;

import java.util.LinkedHashMap;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.bss.Utils;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;

public class AncestralSequenceTrait implements TreeTraitProvider {

	public static String ANCESTRAL_SEQUENCE = "ancestralSequence";
	private Helper helper;

	public AncestralSequenceTrait(
			final LinkedHashMap<NodeRef, int[]> sequenceMap,
			final DataType dataType) {

		helper = new Helper();

		TreeTrait<String> nodeSequence = new TreeTrait.S() {

			@Override
			public String getTraitName() {
				return ANCESTRAL_SEQUENCE;
			}// END: getTraitName

			@Override
			public String getTrait(Tree tree, NodeRef node) {

				int[] intSequence = sequenceMap.get(node);
				
//				System.out.println(intSequence.length);
//				Utils.printArray(intSequence);
				
//				System.exit(-1);
				
				//TODO: null pointer
				String sequence = "AAA";
//				Utils.intArray2Sequence(new Taxon("fake"), //
//						intSequence, //
//						BeagleSequenceSimulator.gapFlag, //
//						dataType).getSequenceString();

				return sequence;
			}// END: getTrait

			@Override
			public Intent getIntent() {
				return Intent.NODE;
			}

		};

		helper.addTrait(nodeSequence);

	}// END: Constructor

	@Override
	public TreeTrait[] getTreeTraits() {

		return helper.getTreeTraits();
	}// END: getTreeTraits

	@Override
	public TreeTrait getTreeTrait(String key) {

		return helper.getTreeTrait(key);
	}// END: getTreeTrait

}// END: class
