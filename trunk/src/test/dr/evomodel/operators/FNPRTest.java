/**
 * 
 */
package test.dr.evomodel.operators;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

import junit.framework.TestCase;

import org.junit.Before;
import org.junit.Test;

import dr.evolution.io.NewickImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.FNPR;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.OperatorFailedException;

/**
 * @author shhn001
 *
 */
public class FNPRTest  extends TestCase{

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
		// if you pick A you can reattach it to 3 new branches
		// if you pick B you can reattach it to 3 new branches
		// if you pick {A,B} you can reattach it to 2 new branches
		// if you pick C you can reattach it to 2 new branches
		// if you pick {A,B,C} you can reattach it to 1 new branch
		// if you pick D you can reattach it to 1 new branch
		// total: 1/12 for every new tree
    	
    	System.out.println("Test 1: Forward");

//        String treeMatch = "(((A,C),D),(B,E));";
        String treeMatch = "(((A,C),D),(E,B));";
        
        int count = 0;
        int reps = 100000;
        
        HashMap<String, Boolean> trees = new HashMap<String, Boolean>();

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5);
                FNPR operator = new FNPR(treeModel, 1);
                operator.doOperation();

                String tree = Tree.Utils.newickNoLengths(treeModel);
//System.out.println(tree);
                if (!trees.containsKey(tree)){
                	trees.put(tree, true);
                }
                if (tree.equals(treeMatch)) {
                    count += 1;
                }

            } catch (OperatorFailedException e) {
                e.printStackTrace();
            }

        }
        System.out.println("Number of trees found:\t" + trees.size());
        Set<String> keys = trees.keySet();
        for (String s : keys){
        	System.out.println(s);
        }
        
        double p_1 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_1);
        System.out.println("Number of expected ratio:\t" + 1.0/12.0);
        assertExpectation(1.0/12.0, p_1, reps);
        
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
