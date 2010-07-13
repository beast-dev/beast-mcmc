package test.dr.evomodel;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.BackboneNodeFilter;
import junit.framework.TestCase;
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

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:2.0);");
        tree = (FlexibleTree) importer.importTree(null);

    }

    public void testMe() {

        TreeTraitProvider.Helper treeTraitProvider = new TreeTraitProvider.Helper();

        TreeTrait dummy = new TreeTrait.D() {
            public String getTraitName() {
                return "test";
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Double getTrait(Tree tree, NodeRef node) {
                return 1.0;
            }
        };

        treeTraitProvider.addTrait(new TreeTrait.FilteredD(dummy, new TreeNodeFilter.ExternalInternalNodeFilter(false, true)));

        Taxa taxonList = new Taxa();

//        treeTraitProvider.addTrait(new TreeTrait.FilteredD("new_name", dummy,
//                new BackboneNodeFilter("name", tree, )));

        System.out.println("Hello there!");

        StringBuffer buffer = new StringBuffer();

                        Tree.Utils.newick(tree, tree.getRoot(), false, Tree.BranchLengthType.LENGTHS_AS_TIME,
                        null, // format
                        null, // branchRates,
                                new TreeTraitProvider[] { treeTraitProvider },
                                null, //idMap, 
                                buffer);

        System.out.println("Tree: " + buffer.toString());
    }


    Tree tree;

}
