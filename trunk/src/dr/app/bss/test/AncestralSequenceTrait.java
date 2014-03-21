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
