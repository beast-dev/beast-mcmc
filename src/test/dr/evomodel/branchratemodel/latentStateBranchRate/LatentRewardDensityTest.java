package test.dr.evomodel.branchratemodel.latentStateBranchRate;

import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import dr.inference.markovjumps.TwoStateSericolaSeriesMarkovReward;
import junit.framework.TestCase;
/**
 * @author JT McCrone
 * This class tests the various classes that estimate the reward for a two state
 * model. The TwoStateSericolaSeriesMarkovReward is used in the SericolaLatentStateBranchModel.
 * TwoStateOccupancyMarkovReward is used in the deprecated LatentStateBranchRateModel
 * We compare the cdf, pdf, and conditional probabilities of these classes where applicable.
 */
public class LatentRewardDensityTest extends TestCase{
    
    private final int dim = 2;
    private final double[] r = new double[]{0.0,1.0}; // rewards
    private final double epsilon = 1e-5;
    private double[] Q(double rate, double bias){
        return new double[]{
            -rate * bias, rate * bias,
            rate * (1.0 - bias), -rate * (1.0 - bias)
        };
    }
    private final double totalTime = 100.0;

    public void testPdfEqualRates(){
        double rate = 0.05;
        double bias = 0.5;
        TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias);
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

        // pdf
        double[] latentTimes =  new double[10];

        for(int i=0; i<latentTimes.length; i++){
            latentTimes[i] = totalTime * i / 10.0;
            double TSOMRPdf = TSOMR.computePdf(latentTimes[i], totalTime, 0, 0);
            double TSSMRPdf = TSSMR.computePdf(latentTimes[i], totalTime, 0, 0);
            double SSMRPdf = SSMR.computePdf(latentTimes[i], totalTime, 0, 0);
            assertEquals(TSOMRPdf, TSSMRPdf,epsilon);
            assertEquals(TSOMRPdf, SSMRPdf,epsilon);
        }


    }
        public void testConditionalEqualRates(){
        double rate = 0.05;
        double bias = 0.5;
        TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias);
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

        //conditional

            double TSOMRConditional = TSOMR.computeConditionalProbability(totalTime, 0, 0);
            double TSSMRConditional = TSSMR.computeConditionalProbability( totalTime, 0, 0);
            double SSMRConditional = SSMR.computeConditionalProbability( totalTime, 0, 0);
            assertEquals(TSOMRConditional, TSSMRConditional,epsilon);
            assertEquals(TSOMRConditional, SSMRConditional,epsilon);

    }
        public void testCdfEqualRates(){
        double rate = 0.05;
        double bias = 0.5;
        // TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias); // no Cdf
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

        double[] latentTimes =  new double[10];
        //conditional
        for(int i=0; i<latentTimes.length; i++){

            double TSSMRCdf = TSSMR.computeCdf(latentTimes[i], totalTime, 0, 0);
            double SSMRCdf = SSMR.computeCdf(latentTimes[i], totalTime, 0, 0);
            assertEquals(SSMRCdf, TSSMRCdf,epsilon);
        }
    }


    public void testPdfUnEqualRates(){
        double rate = 0.05;
        double bias = 0.7;
        TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias);
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

        // pdf
        double[] latentTimes =  new double[10];

        for(int i=0; i<latentTimes.length; i++){
            latentTimes[i] = totalTime * i / 10.0;
            double TSOMRPdf = TSOMR.computePdf(latentTimes[i], totalTime, 0, 0);
            double TSSMRPdf = TSSMR.computePdf(latentTimes[i], totalTime, 0, 0);
            double SSMRPdf = SSMR.computePdf(latentTimes[i], totalTime, 0, 0);
            assertEquals(TSOMRPdf, TSSMRPdf,epsilon);
            assertEquals(TSOMRPdf, SSMRPdf,epsilon);
        }


    }
        public void testConditionalUnEqualRates(){
        double rate = 0.05;
        double bias = 0.7;
        TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias);
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

            double TSOMRConditional = TSOMR.computeConditionalProbability(totalTime, 0, 0);
            double TSSMRConditional = TSSMR.computeConditionalProbability( totalTime, 0, 0);
            double SSMRConditional = SSMR.computeConditionalProbability( totalTime, 0, 0);
            assertEquals(TSOMRConditional, TSSMRConditional,epsilon);
            assertEquals(TSOMRConditional, SSMRConditional,epsilon);

    }
        public void testCdfUnEqualRates(){
        double rate = 0.05;
        double bias = 0.7;
        // TwoStateOccupancyMarkovReward TSOMR = createTwoStateOccupancyMarkovReward(rate, bias); // no Cdf
        TwoStateSericolaSeriesMarkovReward TSSMR = createTwoStateSericolaSeriesMarkovReward(rate, bias);
        SericolaSeriesMarkovReward SSMR = createSericolaSeriesMarkovReward(rate, bias);

        double[] latentTimes =  new double[10];
        //conditional
        for(int i=0; i<latentTimes.length; i++){

            double TSSMRCdf = TSSMR.computeCdf(latentTimes[i], totalTime, 0, 0);
            double SSMRCdf = SSMR.computeCdf(latentTimes[i], totalTime, 0, 0);
            assertEquals(SSMRCdf, TSSMRCdf,epsilon);
        }
    }


    private TwoStateOccupancyMarkovReward createTwoStateOccupancyMarkovReward(double rate , double bias ){
         return new TwoStateOccupancyMarkovReward(Q(rate,bias));
    }
    // this one is used by SericolaLatentBranchModel
    private TwoStateSericolaSeriesMarkovReward createTwoStateSericolaSeriesMarkovReward(double rate , double bias ){
        return new TwoStateSericolaSeriesMarkovReward(Q(rate,bias), r, dim);
    }

    private SericolaSeriesMarkovReward createSericolaSeriesMarkovReward(double rate , double bias ){
        return new  SericolaSeriesMarkovReward(Q(rate,bias), r, dim);
    }
}
