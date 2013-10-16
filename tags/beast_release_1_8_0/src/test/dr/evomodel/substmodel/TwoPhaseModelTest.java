package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.evolution.datatype.Microsatellite;
import dr.evomodel.substmodel.*;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * JUnit test for TwoPhaseModel of microsatellites
 *
 */
public class TwoPhaseModelTest extends TestCase {
    interface Instance {

        public OnePhaseModel getSubModel();

        double getPParam();

        double getMParam();

        double getDistance();

        double[] getPi();

        public double[] getExpectedResult();
    }

    Instance test0 = new Instance() {

        public OnePhaseModel getSubModel(){
            return new AsymmetricQuadraticModel(new Microsatellite(1,6), null);
        }

        public double getPParam(){
            return 0.75;
        }

        public double getMParam(){
            return 0.58;
        }

        public double getDistance(){
            return 0.94;
        }

        public double[] getPi(){
            return new double[]{0.18744637901040673,   0.1607567474109271,   0.15179687357866556,   0.15179687357866556,   0.16075674741092702,    0.18744637901040695};
        }

        public double[] getExpectedResult(){
            return new double[]{
                    0.63398241451779758,    0.23675486586787650,    0.07330431013568975,    0.02919578818459456,    0.01589013386525921,    0.01087248742878196,
                    0.27524752392268953,    0.42098887984631489,    0.19147132450808460,    0.06606140653506744,    0.02837812583210070,    0.01785273935574281,
                    0.09233761983753200,    0.20121070458110243,    0.40820601172807069,    0.19248802732743994,    0.06994493808030816,    0.03581269844554699,
                    0.03581269844554700,    0.06994493808030799,    0.19248802732743994,    0.40820601172807053,    0.20121070458110252,    0.09233761983753204,
                    0.01785273935574278,    0.02837812583210064,    0.06606140653506756,    0.19147132450808424,    0.42098887984631500,    0.27524752392268970,
                    0.01087248742878200,    0.01589013386525917,    0.02919578818459453,    0.07330431013568960,    0.23675486586787683,    0.63398241451779769
            };
        }

    };

    Instance[] all = {test0};


    public void testTwoPhaseModel() {

        for (Instance test : all) {

            OnePhaseModel subModel = test.getSubModel();
            Microsatellite microsat = (Microsatellite)subModel.getDataType();
            Parameter pParam = new Parameter.Default(test.getPParam());
            Parameter mParam = new Parameter.Default(test.getMParam());

            

            TwoPhaseModel tpm = new TwoPhaseModel(microsat, null, subModel, pParam, mParam, null,false);

            int k;

            tpm.computeStationaryDistribution();
            double[] statDist = tpm.getStationaryDistribution();
            final double[] expectedStatDist = test.getPi();
            for (k = 0; k < statDist.length; ++k) {
                assertEquals(statDist[k], expectedStatDist[k], 1e-10);
            }

            int stateCount = microsat.getStateCount();
            double[] mat = new double[stateCount*stateCount];
            tpm.getTransitionProbabilities(test.getDistance(), mat);
            final double[] result = test.getExpectedResult();


            for (k = 0; k < mat.length; ++k) {
                assertEquals(result[k], mat[k], 5e-9);
                //System.out.print(" " + (mat[k]));// - result[k]));
            }

            k = 0;
            for(int i = 0; i < microsat.getStateCount(); i ++){
                for(int j = 0; j < microsat.getStateCount(); j ++){
                    assertEquals(result[k++], tpm.getOneTransitionProbabilityEntry(test.getDistance(), i , j), 1e-10);

                }
            }

            for(int j = 0; j < microsat.getStateCount();j ++){
                double[] colTransitionProb = tpm.getColTransitionProbabilities(test.getDistance(), j);
                for(int i =0 ; i < microsat.getStateCount(); i++){
                    assertEquals(result[i*microsat.getStateCount()+j], colTransitionProb[i], 1e-10);
                }
            }



            for(int i = 0; i < microsat.getStateCount();i ++){
                double[] rowTransitionProb = tpm.getRowTransitionProbabilities(test.getDistance(), i);
                for(int j =0 ; j < microsat.getStateCount(); j++){
                    assertEquals(result[i*microsat.getStateCount()+j], rowTransitionProb[j], 1e-10);
                }
            }
        }

    }


}
