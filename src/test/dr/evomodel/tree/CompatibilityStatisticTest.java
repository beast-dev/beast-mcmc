package test.dr.evomodel.tree;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.CompatibilityStatistic;
import junit.framework.TestCase;

import java.io.IOException;

public class CompatibilityStatisticTest extends TestCase {
    public void setUp() throws Exception {
        super.setUp();

        NewickImporter importer = new NewickImporter("(((1:1.0,2:1.0):0.1,3:1.0):1.0,4:1.0);");
        NewickImporter constraintsImporter = new NewickImporter("((1:1.0,3:1.0,2:1.0):1.0,4:1.0);");

        Tree tree = importer.importTree(null);
        Tree constraintsTree = constraintsImporter.importTree(null);

        compatibilityStatistic = new CompatibilityStatistic("MYSTAT",tree,constraintsTree);
    }

    public void testShouldPass() throws IOException, Importer.ImportException {
        NewickImporter importer2 = new NewickImporter("(((1:1.0,3:1.0):0.1,2:1.0):1.0,4:1.0);");
        Tree tree2 = importer2.importTree(compatibilityStatistic.getTree());
        compatibilityStatistic.setTree(tree2);
        assertTrue(compatibilityStatistic.getBoolean(1));
    }

    public void testShouldFail() throws IOException, Importer.ImportException {
        NewickImporter importer2 = new NewickImporter("(((1:1.0,4:1.0):0.1,3:1.0):1.0,2:1.0);");
        Tree tree2 = importer2.importTree(compatibilityStatistic.getTree());
        compatibilityStatistic.setTree(tree2);
        assertFalse(compatibilityStatistic.getBoolean(1));
    }



    private CompatibilityStatistic compatibilityStatistic;
}
