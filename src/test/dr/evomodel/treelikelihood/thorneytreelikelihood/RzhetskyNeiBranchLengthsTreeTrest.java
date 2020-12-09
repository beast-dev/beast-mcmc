package test.dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.io.NewickImporter;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.treelikelihood.thorneytreelikelihood.BranchLengthProvider;
import dr.evomodel.treelikelihood.thorneytreelikelihood.RzhetskyNeiBranchLengthProvider;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import junit.framework.TestCase;

public class RzhetskyNeiBranchLengthsTreeTrest extends TestCase {

//                          +------------------------------------------------ 0
// +----------------------- |
// |                        +------------------------------------------------ 1
// |
// +-------------------------------------------------------------------------------------------------- 2


    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(1);

        NewickImporter importer = new NewickImporter("((0:1,1:1):0.5,2:2);");
        timeTree= new BigFastTreeModel( importer.importTree(null));

        m = new DistanceMatrix(timeTree);
        m.setElement(0,0,0);
        m.setElement(0,1,2);
        m.setElement(0,2,3.5);

        m.setElement(1,0,2);
        m.setElement(1,1,0);
        m.setElement(1,2,3.5);

        m.setElement(2,0,3.5);
        m.setElement(2,1,3.5);
        m.setElement(2,2,0);


    }

    public void testItWorks(){
        BranchLengthProvider branchLengthProvider = new RzhetskyNeiBranchLengthProvider(m,timeTree);

        System.out.println(branchLengthProvider.getBranchLength(timeTree, timeTree.getNode(0)));

    }

    private TreeModel timeTree;
    private DistanceMatrix m;
}
