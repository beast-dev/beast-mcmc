package test.dr.evolution;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.*;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.BackboneNodeFilter;
import test.dr.inference.trace.TraceCorrelationAssert;

/**
 * @author Marc A. Suchard
 */
public class FilteredTreeTraitTest extends TraceCorrelationAssert {

    public FilteredTreeTraitTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

//        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:2.0);");
//        tree = importer.importTree(null);
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        tree = createPrimateTreeModel();

        treeTraitProvider = new TreeTraitProvider.Helper();

        dummyTrait = new TreeTrait.D() {
            public String getTraitName() {
                return "one";
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Double getTrait(Tree tree, NodeRef node) {
                return 1.0;
            }
        };
    }

    public void testInternalFilter() {

        TreeTrait filteredTrait = new TreeTrait.FilteredD(dummyTrait,
                new TreeNodeFilter.ExternalInternalNodeFilter(false, true));
        treeTraitProvider.addTrait(filteredTrait);
        TreeTrait sumTrait = new TreeTrait.SumOverTreeD(filteredTrait);
        treeTraitProvider.addTrait(sumTrait);

        System.out.println("InternalFilter Test");

        StringBuffer buffer = new StringBuffer();

        TreeUtils.newick(tree, tree.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                null, // format
                null, // branchRates,
                new TreeTraitProvider[]{treeTraitProvider},
                null, //idMap,
                buffer);

        System.out.println("Tree: " + buffer.toString());
        double traitValue = (Double) sumTrait.getTrait(tree, null);
        System.out.println("Trait: " + traitValue);
        assertEquals(traitValue, 5.0);
    }

    public void testBackboneFilter() {

        Taxa taxonList = new Taxa();
        taxonList.addTaxon(taxa[0]);
        taxonList.addTaxon(taxa[1]);

        TreeTrait backboneFilter = new TreeTrait.FilteredD(dummyTrait,
                new BackboneNodeFilter("backbone", tree, taxonList, true, true));
        treeTraitProvider.addTrait(backboneFilter);
        TreeTrait sumTrait = new TreeTrait.SumOverTreeD(backboneFilter);
        treeTraitProvider.addTrait(sumTrait);

        System.out.println("BackboneFilter Test");

        StringBuffer buffer = new StringBuffer();

        TreeUtils.newick(tree, tree.getRoot(), false, TreeUtils.BranchLengthType.LENGTHS_AS_TIME,
                null, // format
                null, // branchRates,
                new TreeTraitProvider[]{treeTraitProvider},
                null, //idMap,
                buffer);

        System.out.println("Tree: " + buffer.toString());
        double traitValue = (Double) sumTrait.getTrait(tree, null);
        System.out.println("Trait: " + traitValue);
        assertEquals(traitValue, 7.0); // TODO Get real result
    }

    Tree tree;
    TreeTrait dummyTrait;
    TreeTraitProvider.Helper treeTraitProvider;

}
