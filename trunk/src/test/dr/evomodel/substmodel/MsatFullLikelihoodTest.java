package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.alignment.Patterns;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.substmodel.AsymmetricQuadraticModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.treelikelihood.TreeLikelihood;

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author Chieh-Hsi Wu
 *
 * JUnit test for testing the full likelihood calculation with microsatellite models.
 */
public class MsatFullLikelihoodTest extends TestCase {

    TreeLikelihood treeLikelihood1;

    public void setUp() throws Exception {
        super.setUp();

        //taxa
        ArrayList<Taxon> taxonList1= new ArrayList<Taxon>();
        Collections.addAll(taxonList1, new Taxon("taxon1"), new Taxon("taxon2"), new Taxon("taxon3"));
        Taxa taxa1 = new Taxa(taxonList1);

        //msat datatype
        Microsatellite msat = new Microsatellite(1,3);
        Patterns msatPatterns = new Patterns(msat, taxa1);
        msatPatterns.addPattern(new int[]{0, 1, 2}); //pattern in the correct code form.

        //create tree
        NewickImporter importer = new NewickImporter("(taxon1:7.5,(taxon2:5.3,taxon3:5.3):2.2);");
        Tree tree =  importer.importTree(null);

        //treeModel
        TreeModel treeModel = new TreeModel(tree);

        //msatsubstModel
        AsymmetricQuadraticModel aqm1 = new AsymmetricQuadraticModel(msat, null);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(aqm1);

        //treeLikelihood
        treeLikelihood1 = new TreeLikelihood(
                    msatPatterns,
                    treeModel,
                    siteModel,
                    null,
                    null,
                    false, false, true, false, false);
    }

    public void testMsatFullLikelihood(){
        double logL1 = -3.29585637705580; //answer calculated by hand.
        assertEquals(logL1, treeLikelihood1.getLogLikelihood(), 1e-10);
        
    }








}
