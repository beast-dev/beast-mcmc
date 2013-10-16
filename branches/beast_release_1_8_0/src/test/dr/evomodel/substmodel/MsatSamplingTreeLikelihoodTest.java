package test.dr.evomodel.substmodel;

import dr.evolution.util.Taxon;
import dr.evolution.util.Taxa;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.alignment.Patterns;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.substmodel.AsymmetricQuadraticModel;
import dr.evomodel.substmodel.TwoPhaseModel;
import dr.evomodel.substmodel.LinearBiasModel;
import dr.evomodel.treelikelihood.MicrosatelliteSamplerTreeLikelihood;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class MsatSamplingTreeLikelihoodTest extends TestCase {

    MicrosatelliteSamplerTreeLikelihood eu1Likelihood;
    MicrosatelliteSamplerTreeLikelihood eu2Likelihood;
    MicrosatelliteSamplerTreeLikelihood ec1Likelihood;
    MicrosatelliteSamplerTreeLikelihood ec2Likelihood;
    MicrosatelliteSamplerTreeLikelihood el1Likelihood;
    MicrosatelliteSamplerTreeLikelihood pu1Likelihood;
    MicrosatelliteSamplerTreeLikelihood pu2Likelihood;
    MicrosatelliteSamplerTreeLikelihood pc1Likelihood;
    MicrosatelliteSamplerTreeLikelihood pc2Likelihood;
    MicrosatelliteSamplerTreeLikelihood pl1Likelihood;


    public void setUp() throws Exception{
        super.setUp();
        //taxa
        ArrayList<Taxon> taxonList3= new ArrayList<Taxon>();
        Collections.addAll(
                taxonList3,
                new Taxon("Taxon1"),
                new Taxon("Taxon2"),
                new Taxon("Taxon3"),
                new Taxon("Taxon4"),
                new Taxon("Taxon5"),
                new Taxon("Taxon6"),
                new Taxon("Taxon7")
        );
        Taxa taxa3 = new Taxa(taxonList3);

        //msat datatype
        Microsatellite msat = new Microsatellite(1,6);
        Patterns msatPatterns = new Patterns(msat, taxa3);
        msatPatterns.addPattern(new int[]{0,1,3,2,4,5,1}); //pattern in the correct code form.

        //create tree
        NewickImporter importer =
                new NewickImporter("(((Taxon1:0.3,Taxon2:0.3):0.6,Taxon3:0.9):0.9,((Taxon4:0.5,Taxon5:0.5):0.3,(Taxon6:0.7,Taxon7:0.7):0.1):1.0);");
        Tree tree =  importer.importTree(null);

        //treeModel
        TreeModel treeModel = new TreeModel(tree);

        //msatsubstModel
        AsymmetricQuadraticModel eu1 = new AsymmetricQuadraticModel(msat, null);


        //create msatSamplerTreeModel
        Parameter internalVal = new Parameter.Default(new double[]{2, 3, 4, 2, 1, 5});
        int[] externalValues = msatPatterns.getPattern(0);
        HashMap<String, Integer> taxaMap = new HashMap<String, Integer>(externalValues.length);
        boolean internalValuesProvided = true;
        for(int i = 0; i < externalValues.length; i++){
            taxaMap.put(msatPatterns.getTaxonId(i),i);
        }
        MicrosatelliteSamplerTreeModel msatTreeModel = new MicrosatelliteSamplerTreeModel("JUnitTestEx", treeModel, internalVal, msatPatterns, externalValues, taxaMap, internalValuesProvided);

        //create msatSamplerTreeLikelihood
        BranchRateModel branchRateModel = new StrictClockBranchRates(new Parameter.Default(1.0));
        eu1Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,eu1, branchRateModel);


        //eu2
        TwoPhaseModel eu2 = new TwoPhaseModel(
                msat,
                null,
                eu1,
                new Parameter.Default(0.0),
                new Parameter.Default(0.4),
                null,
                false
        );
        eu2Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,eu2, branchRateModel);


        //ec1
        LinearBiasModel ec1 = new LinearBiasModel(
                msat,
                null,
                eu1,
                new Parameter.Default(0.48),
                new Parameter.Default(0.0),
                false,
                false,
                false
        );
        ec1Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,ec1, branchRateModel);

        //ec2
        TwoPhaseModel ec2 = new TwoPhaseModel(
                msat,
                null,
                ec1,
                new Parameter.Default(0.0),
                new Parameter.Default(0.4),
                null,
                false
        );
        ec2Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,ec2, branchRateModel);

        //el1
        LinearBiasModel el1 = new LinearBiasModel(
                msat,
                null,
                eu1,
                new Parameter.Default(0.2),
                new Parameter.Default(-0.018),
                true,
                false,
                false
        );
        el1Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,el1, branchRateModel);

        AsymmetricQuadraticModel pu1 = new AsymmetricQuadraticModel(
                msat,
                null,
                new Parameter.Default(1.0),
                new Parameter.Default(0.015),
                new Parameter.Default(0.0),
                new Parameter.Default(1.0),
                new Parameter.Default(0.015),
                new Parameter.Default(0.0),
                false
        );

        pu1Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,pu1, branchRateModel);


        //ec2
        TwoPhaseModel pu2 = new TwoPhaseModel(
                msat,
                null,
                pu1,
                new Parameter.Default(0.0),
                new Parameter.Default(0.4),
                null,
                false
        );
        pu2Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,pu2, branchRateModel);

        //ec1
        LinearBiasModel pc1 = new LinearBiasModel(
                msat,
                null,
                pu1,
                new Parameter.Default(0.48),
                new Parameter.Default(0.0),
                false,
                false,
                false
        );
        pc1Likelihood = new MicrosatelliteSamplerTreeLikelihood(msatTreeModel,pc1, branchRateModel);


    }

    public void testMsatSamplerTreeLikelihood(){
        double eu1LogLik = -40.83979410253295;
        assertEquals(eu1LogLik, eu1Likelihood.getLogLikelihood(), 1e-10);

        double eu2LogLik = -31.06158432691919;
        assertEquals(eu2LogLik, eu2Likelihood.getLogLikelihood(), 1e-10);

        double ec1LogLik = -40.78984094007343;
        assertEquals(ec1LogLik, ec1Likelihood.getLogLikelihood(), 1e-10);

        double ec2LogLik = -31.0303143412092;
        assertEquals(ec2LogLik, ec2Likelihood.getLogLikelihood(), 1e-10);

        double el1LogLik = -40.8979343964233;
        assertEquals(el1LogLik, el1Likelihood.getLogLikelihood(), 1e-10);

        double pu1LogLik = -40.8068080725352;
        assertEquals(pu1LogLik, pu1Likelihood.getLogLikelihood(), 1e-10);

        double pu2LogLik = -31.07280789202575;
        assertEquals(pu2LogLik, pu2Likelihood.getLogLikelihood(), 1e-10);

        double pc1LogLik = -40.7148747667439;
        assertEquals(pc1LogLik, pc1Likelihood.getLogLikelihood(), 1e-10);

    }


}
