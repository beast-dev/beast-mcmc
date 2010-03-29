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
    TreeLikelihood treeLikelihood2;
    TreeLikelihood treeLikelihood3;

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

        setUpExample2();
        setUpExample3();
    }

    private void setUpExample2()throws Exception{
        //taxa
        ArrayList<Taxon> taxonList2= new ArrayList<Taxon>();
        Collections.addAll(
                taxonList2,
                new Taxon("taxon1"),
                new Taxon("taxon2"),
                new Taxon("taxon3"),
                new Taxon("taxon4"),
                new Taxon("taxon5")
        );
        Taxa taxa2 = new Taxa(taxonList2);

        //msat datatype
        Microsatellite msat = new Microsatellite(1,3);
        Patterns msatPatterns = new Patterns(msat, taxa2);
        msatPatterns.addPattern(new int[]{0, 1, 2, 1, 2}); //pattern in the correct code form.

        //create tree
        NewickImporter importer = new NewickImporter("(((taxon1:1.5,taxon2:1.5):1.5,(taxon3:2.1,taxon4:2.1):0.9):0.7,taxon5:3.7);");
        Tree tree =  importer.importTree(null);

        //treeModel
        TreeModel treeModel = new TreeModel(tree);

        //msatsubstModel
        AsymmetricQuadraticModel aqm2 = new AsymmetricQuadraticModel(msat, null);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(aqm2);

        //treeLikelihood
        treeLikelihood2 = new TreeLikelihood(
                    msatPatterns,
                    treeModel,
                    siteModel,
                    null,
                    null,
                    false, false, true, false, false);
    }

    private void setUpExample3() throws Exception{
        //taxa
        ArrayList<Taxon> taxonList3= new ArrayList<Taxon>();
        Collections.addAll(
                taxonList3,
                new Taxon("taxon1"),
                new Taxon("taxon2"),
                new Taxon("taxon3"),
                new Taxon("taxon4"),
                new Taxon("taxon5"),
                new Taxon("taxon6"),
                new Taxon("taxon7")
        );
        Taxa taxa3 = new Taxa(taxonList3);

        //msat datatype
        Microsatellite msat = new Microsatellite(1,4);
        Patterns msatPatterns = new Patterns(msat, taxa3);
        msatPatterns.addPattern(new int[]{0,3,1,2,3,0,1}); //pattern in the correct code form.

        //create tree
        NewickImporter importer =
                new NewickImporter("(((taxon1:0.3,taxon2:0.3):0.6,taxon3:0.9):0.9,((taxon4:0.5,taxon5:0.5):0.3,(taxon6:0.7,taxon7:0.7):0.1):1.0);");
        Tree tree =  importer.importTree(null);

        //treeModel
        TreeModel treeModel = new TreeModel(tree);

        //msatsubstModel
        AsymmetricQuadraticModel aqm3 = new AsymmetricQuadraticModel(msat, null);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(aqm3);

        //treeLikelihood
        treeLikelihood3 = new TreeLikelihood(
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

        double logL2 = -5.51026695214529;
        assertEquals(logL2, treeLikelihood2.getLogLikelihood(), 1e-10);

        double logL3 = -14.39899197302407;
        assertEquals(logL3, treeLikelihood3.getLogLikelihood(), 1e-10);

        
    }








}
