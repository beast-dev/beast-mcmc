/**
 * 
 */
package test.dr.evomodel.operators;

import java.io.IOException;

import dr.evolution.tree.TreeUtils;
import dr.inference.operators.*;
import junit.framework.TestSuite;
import junit.framework.Test;

import dr.evolution.io.NewickImporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.tree.FlexibleTree;
import dr.evomodel.operators.ImportancePruneAndRegraft;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;


/**
 * @author Sebastian Hoehna
 *
 */
public class ImportancePruneAndRegraftTestProblem extends OperatorAssert{

	public static Test suite() {
        return new TestSuite(ImportancePruneAndRegraftTestProblem.class);
    }
	
	/**
	 * Test method for {@link SimpleMCMCOperator#doOperation()}.
	 * @throws ImportException 
	 * @throws IOException 
	 */
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

            TreeModel treeModel = new TreeModel("treeModel", tree5);
            ImportancePruneAndRegraft operator = new ImportancePruneAndRegraft(treeModel, 1.0, 0);
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

            TreeModel treeModel = new TreeModel("treeModel", tree5_2);
            ImportancePruneAndRegraft operator = new ImportancePruneAndRegraft(treeModel, 1.0, 1);
            operator.doOperation();

            String tree = TreeUtils.newickNoLengths(treeModel);

            if (tree.equals(treeMatch)) {
                count += 1;
            }

        }
        double p_2 = (double) count / (double) reps;

        System.out.println("Number of proposals:\t" + count);
        System.out.println("Number of tries:\t" + reps);
        System.out.println("Number of ratio:\t" + p_2);
        System.out.println("Number of expected ratio:\t" + 0.0833333);
        assertExpectation(0.0833333, p_2, reps);
	}
	
	 public OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

		 Parameter rootParameter = treeModel.createNodeHeightsParameter(true, false, false);
	     Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);

	     ImportancePruneAndRegraft operator = new ImportancePruneAndRegraft(treeModel, 1.0, 1);
	     ScaleOperator scaleOperator = new ScaleOperator(rootParameter, 0.75, AdaptationMode.ADAPTATION_ON, 1.0);
	     UniformOperator uniformOperator = new UniformOperator(internalHeights, 1.0);

	     OperatorSchedule schedule = new SimpleOperatorSchedule();
	     schedule.addOperator(operator);
	     schedule.addOperator(scaleOperator);
	     schedule.addOperator(uniformOperator);

	     return schedule;
	 }
}
