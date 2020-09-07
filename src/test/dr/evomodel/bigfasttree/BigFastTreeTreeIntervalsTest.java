package test.dr.evomodel.bigfasttree;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.bigfasttree.CladeAwareSubtreePruneRegraft;
import dr.evomodel.bigfasttree.CladeNodeModel;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.operators.NodeHeightOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import junit.framework.TestCase;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BigFastTreeTreeIntervalsTest extends TestCase {
    public void setUp() throws Exception {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        MathUtils.setSeed(7);
        tree = new BigFastTreeModel(importer.importTree(null));
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

    public void testsmallTree() throws TreeUtils.MissingTaxonException {
        IntervalList bigFastIntervals = new BigFastTreeIntervals(tree);
        bigFastIntervals.calculateIntervals();
        double[] intervals = {0, 0.5, 0.49, 0.01, 0.5, 0, 0.5, 0.5, 0.5, 1.0 };
        //   node             1-2, 4  , 5  , 6  , 0  , 3, 7  , 9  , 8  , 10
        boolean pass = true;
        for (int j = 0; j < bigFastIntervals.getIntervalCount(); j++) {
            if(Math.abs(bigFastIntervals.getInterval(j)- intervals[j])>1E-3){
                System.out.println(bigFastIntervals.getInterval(j)- intervals[j]);
                System.out.println("expected: "+ intervals[j] + " got: "+ bigFastIntervals.getInterval(j));
                pass =false;
                break;
            }
        }
    assertTrue(pass);
    }

    public void testHandelHeightChange() throws TreeUtils.MissingTaxonException {
        IntervalList bigFastIntervals = new BigFastTreeIntervals(tree);
        bigFastIntervals.calculateIntervals();
        double[] intervals = {0, 0.5, 0.49, 0.01, 0.5, 0, 0.5, 0.5, 0.5, 1.5 };
        //   node             1-2, 4  , 5  , 6  , 0  , 3, 7  , 9  , 8  , 10
        boolean pass = true;
        tree.beginTreeEdit();
        tree.setNodeHeight(tree.getNode(10),4.5);
        tree.endTreeEdit();
        bigFastIntervals.calculateIntervals();
        for (int j = 0; j < bigFastIntervals.getIntervalCount(); j++) {
            if(Math.abs(bigFastIntervals.getInterval(j)- intervals[j])>1E-3){
                System.out.println(bigFastIntervals.getInterval(j)- intervals[j]);
                System.out.println("expected: "+ intervals[j] + " got: "+ bigFastIntervals.getInterval(j));
                pass =false;
                break;
            }
        }
        assertTrue(pass);
    }

    //                                                                    0.5
//                                              1.0           +------ 0 * [1.5]
//                                   +----------------------- |(7) *  [2]
//                     1.0           |                        |                                             +-- 1 * [0]
//          +----------------------- |(8) *  [3]              +---------------------------------------------|(6) [0.1]
//          |                        |                                                                      +--- 2 * [0]
//          |                        |                1.5
//          |(10) &  [4]             +-------------------------------- 3 * [1.5]
//          |                                                       2.0
//          |        1.5                   +----------------------------------------------------- 4 & [0.5]
//          +-----------------------------|(9) & [2.5]     1.51
//                                        +-------------------------------------------- 5 & [0.99]
    public void testOrderChange() throws TreeUtils.MissingTaxonException, IOException, Importer.ImportException {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        tree = new BigFastTreeModel(importer.importTree(null));
        IntervalList bigFastIntervals = new BigFastTreeIntervals(tree);
        IntervalList intervals = new TreeIntervals(tree, null, null);
        bigFastIntervals.calculateIntervals();

        //   node             1-2,6,  4  , 5  ,  0  ,  3, 7  , 9  , 8  , 10
        boolean pass = true;
        NodeRef parent = tree.getNode(6);
        NodeRef sibling = tree.getNode(2);
        NodeRef grandParent = tree.getParent(parent);
        NodeRef j = tree.getNode(2);

        tree.beginTreeEdit();

        tree.removeChild(grandParent, parent);
        tree.removeChild(parent, sibling);
        tree.addChild(grandParent, sibling);

        NodeRef jParent = tree.getParent(j);

        tree.removeChild(jParent, j);
        tree.addChild(jParent, parent);
        tree.addChild(parent,j);

        tree.setNodeHeight(parent, 0.1);
        tree.endTreeEdit();

        bigFastIntervals.calculateIntervals();

        List<Integer> missed =new ArrayList<>();
        for (int i = 0; i< bigFastIntervals.getIntervalCount(); i++) {
            if(bigFastIntervals.getInterval(i)!=intervals.getInterval(i)){
                System.out.println("expected: "+ intervals.getInterval(i) + " got: "+ bigFastIntervals.getInterval(i));
                missed.add(i);
            }
        }
        System.out.println(missed);
        assertEquals(0, missed.size());
    }

    public void testCompareIntervals() throws TreeUtils.MissingTaxonException, IOException, Importer.ImportException {
        NewickImporter importer = new NewickImporter("(Lishui/LS557/2020:0,((Netherlands/Utrecht_10015/2020:0.00006795400000000001,USA/IL-NM073/2020:0.00006799599999999999):0.000033976,England/LOND-D604F/2020:0.000101963):0.000033968,Guangdong/2020XN4459-P0041/2020:0.000000005,(Portugal/PT0063/2020:0,(Spain/Zaragoza2486/2020:0.000102605,Scotland/CVR746/2020:0.000000005,Spain/COV000882/2020:0.000067956,Colombia/INS-79253/2020:0.000101944,Uruguay/UY-4/2020:0.000031515):0.000033979,(Spain/CastillaLaMancha201329/2020:0.000000005,Netherlands/NoordHolland_10011/2020:0.000033987):0.00006799,England/LIVE-9CE87/2020:0.00013727299999999998,Spain/Granada-COV002916/2020:0.000033979999999999997):0.000033968,((USA/VI-CDC-3705/2020:0.000000005,Australia/VIC229/2020:0,USA/MA-MGH-00063/2020:0,(USA/WA-S41/2020:0.000068895,USA/WA-UW114/2020:0.000067978,USA/WA-UW17/2020:0.000000005,(USA/WA-S582/2020:0,USA/WA-UW-1682/2020:0.000000005,USA/WA-S994/2020:0.000101934):0.000033955,USA/WA-S121/2020:0.000000005,USA/WA-S154/2020:0.000067982,USA/WA-UW37/2020:0,USA/WA-S321/2020:0,USA/WA-S445/2020:0,USA/WA-S512/2020:0,USA/WA-S33/2020:0.000033979,Canada/BC_6981299/2020:0.000033972,USA/WA-UW-1294/2020:0.000033972,USA/WA-UW-2247/2020:0.000033988,Australia/VIC140/2020:0.000033984,USA/WA-UW61/2020:0.000033972,Canada/BC_8606204/2020:0.000166157,(USA/WA-S734/2020:0,USA/WA-S844/2020:0.000033983):0.000067965,(USA/WA-S1191/2020:0.000067947,USA/WA-S951/2020:0.000101914):0.000095803,Australia/NSW99/2020:0.000101953,(USA/WA-S317/2020:0.000000005,USA/WA-S721/2020:0.00003397):0.00010195700000000001,USA/WA-UW139/2020:0.000135916,USA/WA-S572/2020:0.000033979999999999997,USA/WA-S279/2020:0.000033972,USA/WA-UW28/2020:0.000034002,USA/WA-S114/2020:0.000033969,(USA/WA-S852/2020:0.000203899,(USA/WA-S568/2020:0,USA/WA-S791/2020:0.00006794599999999999):0.000033983):0.000101964,USA/WA-S842/2020:0.000067951):0.000033986,Singapore/302/2020:0.000101947):0.00016677,(((USA/IL-NM0112/2020:0.00003397,USA/IL-NM053/2020:0.000034229,USA/IL-NM059/2020:0.000101967):0.00003397,USA/WI-UW-218/2020:0.000033995):0.000030539,(USA/UT-QDX-63/2020:0,USA/CA-QDX-111/2020:0,USA/TX-HMH0427/2020:0.000203861):0.000101955):0.00023787300000000002):0.000033959,(((Scotland/CVR3203/2020:0.000000005,Scotland/CVR2246/2020:0.000000005,Scotland/GCVR-1714B2/2020:0.000033975999999999995,Scotland/CVR3514/2020:0.000068628):0.000067954,Australia/NT08/2020:0.000034000999999999995):0.00003397,Spain/COV001440/2020:0,Spain/Alcaniz2449/2020:0.000068985,Spain/COV001548/2020:0,USA/WI-WSLH-200057/2020:0.000000005,Spain/Valencia6/2020:0.0000343,Spain/Granada-COV002944/2020:0.000000005,Spain/COV001929/2020:0.000000005,Spain/COV002049/2020:0.000000005,(Spain/Valencia59/2020:0,Spain/Valencia306/2020:0.000000005):0.000033996,(Spain/COV001117/2020:0.00010265,Spain/COV002055/2020:0.000000005,England/20126000104/2020:0.00006758400000000001):0.000033997,Spain/COV001576/2020:0.000000005,Chile/Santiago-1/2020:0.000000005,Spain/COV000721/2020:0.000000005,(Spain/COV001575/2020:0,Spain/COV001505/2020:0):0.000067968,Spain/Madrid_H12_28/2020:0.000067957,Spain/COV001568/2020:0.000033975,England/CAMB-83357/2020:0.000068619,(Spain/Almeria-COV002842/2020:0.000000005,Spain/Malaga-COV002841/2020:0.000000005):0.000067851):0.000169854,Spain/Madrid_LP16_6193/2020:0.00006795299999999999,Singapore/51/2020:0.000044697,(Thailand/Nonthaburi_193/2020:0,Thailand/Bangkok_237/2020:0,Thailand/Bangkok_238/2020:0,((Thailand/Bangkok-0034/2020:0.000000005,Thailand/Bangkok_2295/2020:0,Thailand/Bangkok-0065/2020:0.000047826,Thailand/Bangkok-CONI-0147/2020:0.000033997):0.000033983,Thailand/SI202769-NT/2020:0.000203899):0.000101951):0.00006797500000000001,Shenzhen/SZTH-002/2020:0.000033999);");
        NewickImporter constraintsImporter = new NewickImporter("(Lishui/LS557/2020:0,((Netherlands/Utrecht_10015/2020:0.00006795400000000001,USA/IL-NM073/2020:0.00006799599999999999):0.000033976,England/LOND-D604F/2020:0.000101963):0.000033968,Guangdong/2020XN4459-P0041/2020:0.000000005,(Portugal/PT0063/2020:0,(Spain/Zaragoza2486/2020:0.000102605,Scotland/CVR746/2020:0.000000005,Spain/COV000882/2020:0.000067956,Colombia/INS-79253/2020:0.000101944,Uruguay/UY-4/2020:0.000031515):0.000033979,(Spain/CastillaLaMancha201329/2020:0.000000005,Netherlands/NoordHolland_10011/2020:0.000033987):0.00006799,England/LIVE-9CE87/2020:0.00013727299999999998,Spain/Granada-COV002916/2020:0.000033979999999999997):0.000033968,((USA/VI-CDC-3705/2020:0.000000005,Australia/VIC229/2020:0,USA/MA-MGH-00063/2020:0,(USA/WA-S41/2020:0.000068895,USA/WA-UW114/2020:0.000067978,USA/WA-UW17/2020:0.000000005,(USA/WA-S582/2020:0,USA/WA-UW-1682/2020:0.000000005,USA/WA-S994/2020:0.000101934):0.000033955,USA/WA-S121/2020:0.000000005,USA/WA-S154/2020:0.000067982,USA/WA-UW37/2020:0,USA/WA-S321/2020:0,USA/WA-S445/2020:0,USA/WA-S512/2020:0,USA/WA-S33/2020:0.000033979,Canada/BC_6981299/2020:0.000033972,USA/WA-UW-1294/2020:0.000033972,USA/WA-UW-2247/2020:0.000033988,Australia/VIC140/2020:0.000033984,USA/WA-UW61/2020:0.000033972,Canada/BC_8606204/2020:0.000166157,(USA/WA-S734/2020:0,USA/WA-S844/2020:0.000033983):0.000067965,(USA/WA-S1191/2020:0.000067947,USA/WA-S951/2020:0.000101914):0.000095803,Australia/NSW99/2020:0.000101953,(USA/WA-S317/2020:0.000000005,USA/WA-S721/2020:0.00003397):0.00010195700000000001,USA/WA-UW139/2020:0.000135916,USA/WA-S572/2020:0.000033979999999999997,USA/WA-S279/2020:0.000033972,USA/WA-UW28/2020:0.000034002,USA/WA-S114/2020:0.000033969,(USA/WA-S852/2020:0.000203899,(USA/WA-S568/2020:0,USA/WA-S791/2020:0.00006794599999999999):0.000033983):0.000101964,USA/WA-S842/2020:0.000067951):0.000033986,Singapore/302/2020:0.000101947):0.00016677,(((USA/IL-NM0112/2020:0.00003397,USA/IL-NM053/2020:0.000034229,USA/IL-NM059/2020:0.000101967):0.00003397,USA/WI-UW-218/2020:0.000033995):0.000030539,(USA/UT-QDX-63/2020:0,USA/CA-QDX-111/2020:0,USA/TX-HMH0427/2020:0.000203861):0.000101955):0.00023787300000000002):0.000033959,(((Scotland/CVR3203/2020:0.000000005,Scotland/CVR2246/2020:0.000000005,Scotland/GCVR-1714B2/2020:0.000033975999999999995,Scotland/CVR3514/2020:0.000068628):0.000067954,Australia/NT08/2020:0.000034000999999999995):0.00003397,Spain/COV001440/2020:0,Spain/Alcaniz2449/2020:0.000068985,Spain/COV001548/2020:0,USA/WI-WSLH-200057/2020:0.000000005,Spain/Valencia6/2020:0.0000343,Spain/Granada-COV002944/2020:0.000000005,Spain/COV001929/2020:0.000000005,Spain/COV002049/2020:0.000000005,(Spain/Valencia59/2020:0,Spain/Valencia306/2020:0.000000005):0.000033996,(Spain/COV001117/2020:0.00010265,Spain/COV002055/2020:0.000000005,England/20126000104/2020:0.00006758400000000001):0.000033997,Spain/COV001576/2020:0.000000005,Chile/Santiago-1/2020:0.000000005,Spain/COV000721/2020:0.000000005,(Spain/COV001575/2020:0,Spain/COV001505/2020:0):0.000067968,Spain/Madrid_H12_28/2020:0.000067957,Spain/COV001568/2020:0.000033975,England/CAMB-83357/2020:0.000068619,(Spain/Almeria-COV002842/2020:0.000000005,Spain/Malaga-COV002841/2020:0.000000005):0.000067851):0.000169854,Spain/Madrid_LP16_6193/2020:0.00006795299999999999,Singapore/51/2020:0.000044697,(Thailand/Nonthaburi_193/2020:0,Thailand/Bangkok_237/2020:0,Thailand/Bangkok_238/2020:0,((Thailand/Bangkok-0034/2020:0.000000005,Thailand/Bangkok_2295/2020:0,Thailand/Bangkok-0065/2020:0.000047826,Thailand/Bangkok-CONI-0147/2020:0.000033997):0.000033983,Thailand/SI202769-NT/2020:0.000203899):0.000101951):0.00006797500000000001,Shenzhen/SZTH-002/2020:0.000033999);");

//        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
//        NewickImporter constraintsImporter = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");

        tree = new BigFastTreeModel(importer.importTree(null));
        Tree constraints = constraintsImporter.importTree(null);

        CladeNodeModel clades = new CladeNodeModel(constraints, tree);

        IntervalList intervals = new TreeIntervals(tree,null,null);
        BigFastTreeIntervals bigFastTreeIntervals = new BigFastTreeIntervals(tree);

        CladeAwareSubtreePruneRegraft op = new CladeAwareSubtreePruneRegraft(clades, 1,1);
        NodeHeightOperator nh = new NodeHeightOperator(tree,1,1000000, NodeHeightOperator.OperatorType.UNIFORM, AdaptationMode.ADAPTATION_OFF,0.25);
        NodeHeightOperator root = new NodeHeightOperator(tree,1,0.75, NodeHeightOperator.OperatorType.SCALEROOT, AdaptationMode.ADAPTATION_OFF,0.25);
        boolean pass = true;

        MathUtils.setSeed(2);
        for (int i = 0; i < 10000; i++) {
            if(MathUtils.nextBoolean()){
//                op.doOperation();
                root.doOperation();
            }else{
                nh.doOperation();
            }
            intervals.calculateIntervals();
//            bigFastIntervals.makeDirty();
            bigFastTreeIntervals.calculateIntervals();
            for (int j = 0; j < bigFastTreeIntervals.getIntervalCount(); j++) {
                if(intervals.getInterval(j)!= bigFastTreeIntervals.getInterval(j)){
                    System.out.println(i);
                    System.out.println("interval wrong");
                    pass =false;
                    break;
                }
            }
            for (int j = 0; j < bigFastTreeIntervals.getIntervalCount(); j++) {
                if(intervals.getLineageCount(j)!= bigFastTreeIntervals.getLineageCount(j)){
                    System.out.println(i);
                    System.out.println("lineage Counts wrong: " +j);
                    System.out.println("expected: " + intervals.getLineageCount(j));
                    System.out.println("got " +bigFastTreeIntervals.getLineageCount(j) );

                    pass =false;
                    break;
                }
            }
            for (int j = 0; j < bigFastTreeIntervals.getIntervalCount(); j++) {
                if(intervals.getIntervalTime(j)!= bigFastTreeIntervals.getIntervalTime(j)){
                    System.out.println(i);
                    System.out.println("times wrong");
                    pass =false;
                    break;
                }
            }
            if (!pass) {
                break;
            }
        }
        assertTrue(pass);
    }


    private TreeModel tree;
}
