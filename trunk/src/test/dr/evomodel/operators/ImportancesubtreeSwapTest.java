/**
 * 
 */
package test.dr.evomodel.operators;


import java.io.IOException;

import junit.framework.TestCase;

import org.junit.Before;
import org.junit.Test;

import dr.evolution.io.NewickImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.ImportanceSubtreeSwap;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;

/**
 * @author shhn001
 *
 */
public class ImportancesubtreeSwapTest extends TestCase{


	private FlexibleTree tree5;
    private FlexibleTree tree6;
    
	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		super.setUp();

        NewickImporter importer = new NewickImporter(
                "((((A:1.0,B:1.0):1.0,C:2.0):1.0,D:4.0):1.0,E:5.0);");
        tree5 = (FlexibleTree) importer.importTree(null);

        importer = new NewickImporter(
                "(((((A:1.0,B:1.0):1.0,C:2.0):1.0,D:4.0):1.0,E:5.0),F:6.0);");
        tree6 = (FlexibleTree) importer.importTree(null);
	}

	/**
	 * Test method for {@link dr.evomodel.operators.ImportanceSubtreeSwap#doOperation()}.
	 * @throws ImportException 
	 * @throws IOException 
	 */
	@Test
	public void testDoOperation() throws IOException, ImportException {
		// probability of picking (A,B) node is 1/(2n-3) = 1/7
        // probability of swapping with D is 1/2
        // total = 1/14

        //probability of picking {D} node is 1/(2n-3) = 1/7
        //probability of picking {A,B} is 1/5
        // total = 1/35

        //total = 1/14 + 1/35 = 7/70 = 0.1
    	
    	System.out.println("Test 1: Forward");

        String treeMatch = "(((D,C),(A,B)),E);";
        
        int count = 0;
        int reps = 1000000;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5);
                ImportanceSubtreeSwap operator = new ImportanceSubtreeSwap(treeModel, 1.0, 1);
                operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        double p_1 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_1);
        System.out.println("Number of expected ratio:\t" + 0.1);
        assertExpectation(0.1, p_1, reps);
        
        // lets see what the backward probability is for the hastings ratio
        
        // (((D:2.0,C:2.0):1.0,(A:1.0,B:1.0):2.0):1.0,E:4.0) -> ((((A,B),C),D),E)
        
        // probability of picking (A,B) node is 1/(2n-3) = 1/7
        // probability of swapping with D is 1/3
        // total = 1/21

        //probability of picking {D} node is 1/(2n-2) = 1/7
        //probability of picking {A,B} is 1/4
        // total = 1/28

        //total = 1/21 + 1/28 = 7/84 = 0.08333333
        
    	System.out.println("Test 2: Backward");
        
        treeMatch = "((((A,B),C),D),E);";
        NewickImporter importer = new NewickImporter("(((D:2.0,C:2.0):1.0,(A:1.0,B:1.0):2.0):1.0,E:4.0);");
        FlexibleTree tree5_2 = (FlexibleTree) importer.importTree(null);

        count = 0;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5_2);
                ImportanceSubtreeSwap operator = new ImportanceSubtreeSwap(treeModel, 1.0, 1);
                operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        double p_2 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_2);
        System.out.println("Number of expected ratio:\t" + 0.0833333);
        assertExpectation(0.0833333, p_2, reps);
	}
	
	/**
     * @param ep    the expected (binomial) probability of success
     * @param ap    the actual proportion of successes
     * @param count the number of attempts
     */
    protected void assertExpectation(double ep, double ap, int count) {

        if (count * ap < 5 || count * (1 - ap) < 5) throw new IllegalArgumentException();

        double stdev = Math.sqrt(ap * (1.0 - ap) * count) / count;
        double upper = ap + 2 * stdev;
        double lower = ap - 2 * stdev;

        assertTrue("Expected p=" + ep + " but got " + ap + " +/- " + stdev,
                upper > ep && lower < ep);

    }

}
