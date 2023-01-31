/**
 * 
 */
package test.dr.evomodel.operators;

import java.io.IOException;
import java.util.HashMap;
import java.util.Set;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.operators.*;
import junit.framework.TestSuite;
import junit.framework.Test;

import dr.evolution.io.Importer.ImportException;
import dr.evomodel.operators.FNPR;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

/**
 * @author shhn001
 *
 */
public class FNPRTest extends OperatorAssert{

    public static Test suite() {
        return new TestSuite(FNPRTest.class);
    }
    
	/**
	 * Test method for {@link SimpleMCMCOperator#doOperation()}.
	 * @throws ImportException 
	 * @throws IOException 
	 */
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

            DefaultTreeModel treeModel = new DefaultTreeModel("treeModel", tree5);
            FNPR operator = new FNPR(treeModel, 1);
            operator.doOperation();

            String tree = TreeUtils.newickNoLengths(treeModel);
//System.out.println(tree);
            if (!trees.containsKey(tree)){
                trees.put(tree, true);
            }
            if (tree.equals(treeMatch)) {
                count += 1;
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
	
	public OperatorSchedule getOperatorSchedule(DefaultTreeModel treeModel) {

        Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
        Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);

        FNPR operator = new FNPR(treeModel, 1.0);
        ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, AdaptationMode.ADAPTATION_ON, 1.0);
        UniformOperator uniformOperator = new UniformOperator(internalHeights, 1.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(operator);
        schedule.addOperator(scaleOperator);
        schedule.addOperator(uniformOperator);

        return schedule;
    }

}
