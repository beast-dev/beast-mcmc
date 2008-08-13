package test.dr.evomodel.operators;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * @author Alexei Drummond
 */
public class ExchangeOperatorTest extends TestCase {

    FlexibleTree tree;

    public void setUp() throws Exception {
        super.setUp();

        NewickImporter importer = new NewickImporter(
                "((((A:1.0,B:1.0):1.0,C:2.0),D:3.0):1.0, E:4.0);");
        tree = (FlexibleTree) importer.importTree(null);
    }


    public void testWideExchangeOperator() {

        // probability of picking (A,B) node is 1/(2n-2) = 1/8
        // probability of swapping with D is 1/2
        // total = 1/16

        //probability of picking {D} node is 1/(2n-2) = 1/8
        //probability of picking {A,B} is 1/5
        // total = 1/40

        //total = 1/16 + 1/40 = 0.0625 + 0.025 = 0.085


        String treeMatch = "(((D,C),(A,B)),E);";

        int count = 0;
        int reps = 10000;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree);
                ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
                double logq = operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        double p = (double) count / (double) reps;

        assertExpectation(0.085, p, reps);
    }

    /**
     */
    protected void assertExpectation(double ep, double ap, int count) {

        if (count * ap < 5 || count * (1 - ap) < 5) throw new IllegalArgumentException();

        double stdev = Math.sqrt(ap * (1.0 - ap) * count) / count;
        double upper = ap + 2 * stdev;
        double lower = ap - 2 * stdev;

        assertTrue("Expected p=" + ep + " but got " + ap + " +/- " + stdev,
                upper > ep && lower < ep);

    }


    public static Test suite() {
        return new TestSuite(ExchangeOperatorTest.class);
    }

}
