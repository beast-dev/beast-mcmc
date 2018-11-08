/**
 * 
 */
package test.dr.evomodel.operators;

//import static org.junit.Assert.*;
import java.io.IOException;

import dr.evolution.tree.TreeUtils;
import junit.framework.Test;
import junit.framework.TestSuite;

import dr.evolution.io.Importer.ImportException;
import dr.evomodel.operators.ExchangeOperator;
//import dr.evomodel.operators.NNI;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;

/**
 * @author shhn001
 *
 */
public class NarrowExchangeTest  extends OperatorAssert{
    
	public static Test suite() {
        return new TestSuite(NarrowExchangeTest.class);
    }


	/**
	 * Test method for {@link dr.evomodel.operators.ExchangeOperator#narrow()}.
	 */
	public void testNarrow() throws IOException, ImportException {
		// probability of picking B node is 1/(2n-4) = 1/6
        // probability of swapping it with C is 1/1
        // total = 1/6
    	
    	System.out.println("Test 1: Forward");

        String treeMatch = "((((A,C),B),D),E);";
        
        int count = 0;
        int reps = 100000;

        for (int i = 0; i < reps; i++) {

            TreeModel treeModel = new TreeModel("treeModel", tree5);
            ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 1);
            operator.doOperation();

            String tree = TreeUtils.newickNoLengths(treeModel);

            if (tree.equals(treeMatch)) {
                count += 1;
            }

        }
        double p_1 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_1);
        System.out.println("Number of expected ratio:\t" + 1.0/6.0);
        assertExpectation(1.0/6.0, p_1, reps);
        
	}
	
	public OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);

        ExchangeOperator operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 1.0);
        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, AdaptationMode.ADAPTATION_ON, 1.0);
        UniformOperator uniformOperator = new UniformOperator(internalHeights, 1.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(operator);
        schedule.addOperator(scaleOperator);
        schedule.addOperator(uniformOperator);

        return schedule;
    }
}
