/**
 * 
 */
package test.dr.evomodel.operators;

import java.io.IOException;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
//import dr.evolution.io.NewickImporter;
import dr.evolution.io.Importer.ImportException;
//import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.GibbsSubtreeSwap;
//import dr.evomodel.operators.ImportanceSubtreeSwap;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;


/**
 * @author Sebastian Hoehna
 *
 */
public class GibbsSubtreeSwapTestProblem extends OperatorAssert{


	public static Test suite() {
        return new TestSuite(GibbsSubtreeSwapTestProblem.class);
    }
	
	/**
	 * Test method for {@link dr.evomodel.operators.GibbsSubtreeSwap#doOperation()}.
	 * @throws ImportException 
	 * @throws IOException 
	 */
	public void testDoOperation() throws IOException, ImportException {
		// assumes that the posterior is equal for all trees!!!		
		
		// probability of picking (A,B) node is 1/(2n-3) = 1/7
        // probability of swapping with D is 1/2
        // total = 1/14

        //probability of picking {D} node is 1/(2n-3) = 1/7
        //probability of picking {A,B} is 1/5
        // total = 1/35

        //total = 1/14 + 1/35 = 7/70 = 0.1
		
		// now we calculate the same for the backward proposal
		// this is needed for the Hastings ratio
		
		// probability of picking (A,B) node is 1/(2n-3) = 1/7
        // probability of swapping with D is 1/3
        // total = 1/21

        //probability of picking {D} node is 1/(2n-3) = 1/7
        //probability of picking {A,B} is 1/4
        // total = 1/28

        //total = 1/21 + 1/28 = 4/84 + 3/84 = 7/84 = 1/12 
    	
    	System.out.println("Test 1: Forward");

        String treeMatch = "(((D,C),(A,B)),E);";
        
        int count = 0;
        int reps = 100000;

        for (int i = 0; i < reps; i++) {

            try {
                TreeModel treeModel = new TreeModel("treeModel", tree5);
                GibbsSubtreeSwap operator = new GibbsSubtreeSwap(treeModel, false, 1.0);
                double hr = operator.operate(null, null);

                String tree = Tree.Utils.newickNoLengths(treeModel);

                if (tree.equals(treeMatch)) {
//                	System.out.println("Expected Hastings ratio = " + 5.0/6.0 + " in log = " + Math.log(5.0/6.0));
//                	System.out.println("Obtained Hastings ratio = " + Math.exp(hr) + " in log = " + hr);
                	TestCase.assertEquals(Math.log(5.0/6.0), hr, 0.00000001);
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
	}
	
	 public OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

	        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
	        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);

	        GibbsSubtreeSwap operator = new GibbsSubtreeSwap(treeModel, false, 1.0);
	        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, CoercionMode.COERCION_ON, 1.0);
	        UniformOperator uniformOperator = new UniformOperator(internalHeights, 1.0);

	        OperatorSchedule schedule = new SimpleOperatorSchedule();
	        schedule.addOperator(operator);
	        schedule.addOperator(scaleOperator);
	        schedule.addOperator(uniformOperator);

	        return schedule;
	    }
}
