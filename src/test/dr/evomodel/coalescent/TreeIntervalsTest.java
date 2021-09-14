package test.dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.io.NewickImporter;

import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.Arrays;

public class TreeIntervalsTest extends TestCase {

    private TreeIntervalList treeIntervals;

    public void setUp() throws Exception {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        MathUtils.setSeed(7);
        tree = new DefaultTreeModel(importer.importTree(null));
        treeIntervals = new TreeIntervals(tree,true);
        treeIntervals.calculateIntervals();
    }


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
    public void testsmallTree() {

        double[] intervals = {0, 0.5, 0.49, 0.01, 0.5, 0, 0.5, 0.5, 0.5, 1.0};
        //   node             1-2, 4  , 5  , 6  , 0  , 3, 7  , 9  , 8  , 10
        boolean pass = true;
        for (int j = 0; j < treeIntervals.getIntervalCount(); j++) {
            if (Math.abs(treeIntervals.getInterval(j) - intervals[j]) > 1E-3) {
                System.out.println(treeIntervals.getInterval(j) - intervals[j]);
                System.out.println("expected: " + intervals[j] + " got: " + treeIntervals.getInterval(j));
                pass = false;
                break;
            }
        }
        assertTrue(pass);
    }
    public void testCoalescentNode(){
        int[] coalIntervals = {3,6,7,8,9};
        int[] nodeNumbers  ={6,7,9,8,10};
        boolean pass = true;
        for (int i=0;i<coalIntervals.length;i++) {
            if(treeIntervals.getCoalescentNode(coalIntervals[i]).getNumber()!=nodeNumbers[i]){
                System.out.print(coalIntervals[i]);
                System.out.print("!=");
                System.out.println(nodeNumbers[i]);
                pass=false;
                break;
            }
        }
        assertTrue(pass);
    }

    public void testCoalscentIntervals(){
        double[] trueIntervals = {0.01,0.5,0.5,0.5,1.0};
        double[] calcIntervals = treeIntervals.getCoalescentIntervals();
        boolean pass = true;
        double epsilon = 1e-9;

        assertEquals(trueIntervals.length,calcIntervals.length);
        for (int i=0;i<trueIntervals.length;i++) {
            if (calcIntervals[i] < trueIntervals[i]-epsilon || calcIntervals[i]>trueIntervals[i]+epsilon) {
                System.out.println(Arrays.toString(trueIntervals));
                System.out.println(Arrays.toString(calcIntervals));
                pass=false;
                break;
            }
        }
        assertTrue(pass);
    }
    public void testNodeNumbersForInterval() {
        int[][] intervalNodes= {{1,2},{2, 4} ,{4 , 5},{5, 6},{6  , 0}  ,{0, 3},{3, 7} ,{7 , 9},{9  , 8},{8 , 10}};
        boolean pass = true;
        for (int j = 0; j < treeIntervals.getIntervalCount(); j++) {
            int[] treeIntervalNodes = treeIntervals.getNodeNumbersForInterval(j);
            if (intervalNodes[j][0] != treeIntervalNodes[0] || intervalNodes[j][1] != treeIntervalNodes[1]) {
                System.out.println(Arrays.toString(intervalNodes));
                System.out.println(Arrays.toString(treeIntervalNodes));
                pass = false;
                break;
            }
        }
        assertTrue(pass);
    }

    public void testIntervalsForNode(){

        int[][] intervalNodeIntervals= {{4,5},{0} ,{0 , 1},{5, 6},{1 , 2} ,{2, 3},{3, 4} ,{6 , 7},{8, 9},{7, 8},{9}};
        boolean pass =true;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            int[] treeIntervalNodeIntervals = treeIntervals.getIntervalsForNode(i);
            for (int j = 0; j < intervalNodeIntervals[i].length; j++) {
                if (intervalNodeIntervals[i][j] != treeIntervalNodeIntervals[j]) {
                    System.out.println(Arrays.toString(intervalNodeIntervals[i]));
                    System.out.println(Arrays.toString(treeIntervalNodeIntervals));
                    pass = false;
                    break;
                }
            }
        }
        assertTrue(pass);
    }

    /**
     * This test the internal node order sorting. It takes an array sorted by interval and
     * returns the same array sorted by node number where the node is the node that starts the
     * interval.
     */
    public void testSortByNodeNumber(){
        //                    these are also the node indices for internal nodes
        double[] intervalOrder = {6, 7, 9, 8, 10};
        double[] nodeOrder = {6,7,8,9,10};

        double[] sorted = treeIntervals.sortByNodeNumbers(intervalOrder);
        double epsilon = 1e-9;
        boolean pass=true;
        for (int i = 1; i < sorted.length; i++) {
            if (sorted[i] < nodeOrder[i]-epsilon || sorted[i]>nodeOrder[i]+epsilon) {
                System.out.println(Arrays.toString(sorted));
                System.out.println(Arrays.toString(nodeOrder));
                pass = false;
                break;
            }
        }
        assertTrue(pass);
    }
    public void testHandelHeightChange(){
        // index              0  1     2    3    4     5   6    7   8    9
        double[] intervals = {0, 0.5, 0.49, 0.01, 0.5, 0, 0.5, 0.5, 0.5, 1.5};
        //   node             1-2, 4  , 5  , 6  , 0  , 3, 7  , 9  , 8  , 10
        boolean pass = true;
        tree.beginTreeEdit();
        tree.setNodeHeight(tree.getNode(10), 4.5);
        tree.endTreeEdit();
        treeIntervals.calculateIntervals();
        for (int j = 0; j < treeIntervals.getIntervalCount(); j++) {
            if (Math.abs(treeIntervals.getInterval(j) - intervals[j]) > 1E-3) {
                System.out.println(treeIntervals.getInterval(j) - intervals[j]);
                System.out.println("expected: " + intervals[j] + " got: " + treeIntervals.getInterval(j));
                pass = false;
                break;
            }
        }
        assertTrue(pass);
    }
    private TreeModel tree;

}
