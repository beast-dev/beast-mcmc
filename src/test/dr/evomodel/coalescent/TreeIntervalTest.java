package test.dr.evomodel.coalescent;

import static org.junit.Assert.assertTrue;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.io.NewickImporter;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import junit.framework.TestCase;


public class TreeIntervalTest extends TestCase {

    private dr.evolution.coalescent.TreeIntervals coalescentTreeIntervals;
    private dr.evomodel.coalescent.TreeIntervals modelFastTreeIntervals;
    private dr.evomodel.coalescent.TreeIntervals modelTreeIntervals;
    private dr.evomodel.bigfasttree.BigFastTreeIntervals bigFastTreeIntervals;

    private TreeModel tree;


//    clade * 0-3 root 8
//    clade & 4-5 root 10
//                                                                    0.5
//                                              1.0           +------ 0 * [1.5]
//                                   +----------------------- |(7) *  [2]                          1.0
//                     1.0           |                        |        1.0             +----------------------- 1 * [0]
//          +----------------------- |(8) *  [3]              +----------------------- |(6) [1]    1.0
//          |                        |                                                 +----------------------- 2 * [0]
//          |                        |                1.5
//          |(10) &  [4]             +-------------------------------- 3 * [1.5]
//          |                                                       2.0
//          |        1.5                   +----------------------------------------------------- 4 & [0.5]
//          +-----------------------------|(9) & [2.5]     1.51
//                                        +-------------------------------------------- 5 & [0.99]



    public void setUp() throws Exception {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        MathUtils.setSeed(7);
        tree = new DefaultTreeModel(importer.importTree(null));
        
        bigFastTreeIntervals = new BigFastTreeIntervals(tree);
        bigFastTreeIntervals.calculateIntervals();

        modelFastTreeIntervals = new dr.evomodel.coalescent.TreeIntervals(tree);
        modelFastTreeIntervals.calculateIntervals();

        //building Node mapping uses old treeInvervals not fastIntervals
        modelTreeIntervals = new dr.evomodel.coalescent.TreeIntervals(tree,true);
        modelTreeIntervals.calculateIntervals();

        coalescentTreeIntervals = new dr.evolution.coalescent.TreeIntervals(tree);
        coalescentTreeIntervals.calculateIntervals();
    }
    //Aggregates
    public void testIntervalsCount() {
        assertTrue(bigFastTreeIntervals.getIntervalCount() == modelFastTreeIntervals.getIntervalCount());
        assertTrue(bigFastTreeIntervals.getIntervalCount() == modelTreeIntervals.getIntervalCount());
        assertTrue(bigFastTreeIntervals.getIntervalCount() == coalescentTreeIntervals.getIntervalCount());
    
    }

    public void testGetSampleCount(){
        assertTrue(bigFastTreeIntervals.getSampleCount() == modelFastTreeIntervals.getSampleCount());
        assertTrue(bigFastTreeIntervals.getSampleCount() == modelTreeIntervals.getSampleCount());
        assertTrue(bigFastTreeIntervals.getSampleCount() == coalescentTreeIntervals.getSampleCount());
    }

    public void testGetStartTime(){
        assertTrue(bigFastTreeIntervals.getStartTime() == modelFastTreeIntervals.getStartTime());
        assertTrue(bigFastTreeIntervals.getStartTime() == modelTreeIntervals.getStartTime());
        assertTrue(bigFastTreeIntervals.getStartTime() == coalescentTreeIntervals.getStartTime());
    }

    public void testGetTotalDuration(){
        assertTrue(bigFastTreeIntervals.getTotalDuration() == modelFastTreeIntervals.getTotalDuration());
        assertTrue(bigFastTreeIntervals.getTotalDuration() == modelTreeIntervals.getTotalDuration());
        assertTrue(bigFastTreeIntervals.getTotalDuration() == coalescentTreeIntervals.getTotalDuration());
    }

    public void testIsBinaryCoalescent(){
        assertTrue(bigFastTreeIntervals.isBinaryCoalescent() == modelFastTreeIntervals.isBinaryCoalescent());
        assertTrue(bigFastTreeIntervals.isBinaryCoalescent() == modelTreeIntervals.isBinaryCoalescent());
        assertTrue(bigFastTreeIntervals.isBinaryCoalescent() == coalescentTreeIntervals.isBinaryCoalescent());
    }

    public void testIsCoalescentOnly(){
        assertTrue(bigFastTreeIntervals.isCoalescentOnly() == modelFastTreeIntervals.isCoalescentOnly());
        assertTrue(bigFastTreeIntervals.isCoalescentOnly() == modelTreeIntervals.isCoalescentOnly());
        assertTrue(bigFastTreeIntervals.isCoalescentOnly() == coalescentTreeIntervals.isCoalescentOnly());
    }
}