package dr.app.tools;

import java.io.IOException;

import dr.app.util.Arguments;
import dr.inference.markovjumps.SericolaSeriesMarkovReward;
import dr.inference.markovjumps.TwoStateOccupancyMarkovReward;
import dr.inference.markovjumps.TwoStateSericolaSeriesMarkovReward;

/**
 * The goal is to calculate the pdf and cdf of time spent in a given state 
 * for a continuous-time Markov chain with rewards.
 * 
 * There are several implementations that can achieve this. 
 * We will compare them here and output the results in a json 
 * file format.
 * @author JT McCrone
 */
public class MarkovRewardDensityCalculator {
    private double[] Q ;
    private final int dim = 2;
    private final double[] r = new double[]{0.0,1.0}; // rewards

    private TwoStateOccupancyMarkovReward twoStateOccupancyMarkovReward;
    private TwoStateSericolaSeriesMarkovReward twoStateSericolaSeriesMarkovReward;
    private SericolaSeriesMarkovReward sericolaSeriesMarkovReward;

    private MarkovRewardDensityCalculator(double rate,double bias){
        this.Q = new double[]{
            -rate * bias, rate * bias,
            rate * (1.0 - bias), -rate * (1.0 - bias)
    };

    this.twoStateOccupancyMarkovReward = new TwoStateOccupancyMarkovReward(Q);
    this.twoStateSericolaSeriesMarkovReward = new TwoStateSericolaSeriesMarkovReward(Q, r, dim);
    this.sericolaSeriesMarkovReward = new SericolaSeriesMarkovReward(Q, r, dim);

    }


    void run(double time){
        double s=0.0;
        double[] times =  new double[100];
       
        double pdf;
        double cdf;
        double conditional;
        
        System.out.print("{\"data\":[");

        for(int i=0; i<times.length; i++){
            times[i] = time * i / 100.0;

            System.out.print("{");

//TODO remove duplication
                    pdf = twoStateOccupancyMarkovReward.computePdf(times[i], time, 0, 0);
                    // cdf = twoStateOccupancyMarkovReward.computeCdf(times[i], time, 0, 0);
                    conditional = twoStateOccupancyMarkovReward.computeConditionalProbability( time, 0, 0);
                    
                    System.out.print("\"time\": " + times[i] + ",");
                    System.out.print("\"totalTime\": " + time + ",");
                    // System.out.print("\"twoStateOccupancyMarkovReward\": {");
                    System.out.print("\"pdf\": " + pdf + ",");
                    // System.out.print("\"cdf\": " + cdf + ",");
                    System.out.print("\"conditional\": " + conditional + ",");
                    System.out.print("\"implementation\": \"twoStateOccupancyMarkovReward\"");
                
            System.out.print("},");

            System.out.print("{");
                
                    pdf = twoStateSericolaSeriesMarkovReward.computePdf(times[i], time, 0, 0);
                    cdf = twoStateSericolaSeriesMarkovReward.computeCdf(times[i], time, 0, 0);
                    conditional = twoStateSericolaSeriesMarkovReward.computeConditionalProbability(time, 0, 0);
                    System.out.print("\"time\": " + times[i] + ",");
                    System.out.print("\"totalTime\": " + time + ",");
                    System.out.print("\"pdf\": " + pdf + ",");
                    System.out.print("\"cdf\": " + cdf + ",");
                    System.out.print("\"conditional\": " + conditional + ",");
                    System.out.print("\"implementation\": \"twoStateSericolaSeriesMarkovReward\"");

            System.out.print("},");

            System.out.print("{");
                    pdf = sericolaSeriesMarkovReward.computePdf(times[i], time, 0, 0);
                    cdf = sericolaSeriesMarkovReward.computeCdf(times[i], time, 0, 0);
                    conditional = sericolaSeriesMarkovReward.computeConditionalProbability(time, 0, 0);

                    System.out.print("\"time\": " + times[i] + ",");
                    System.out.print("\"totalTime\": " + time + ",");
                    System.out.print("\"pdf\": " + pdf + ",");
                    System.out.print("\"cdf\": " + cdf + ",");
                    System.out.print("\"conditional\": " + conditional + ",");
                    System.out.print("\"implementation\": \"sericolaSeriesMarkovReward\"");

                if(i < times.length - 1){
                    System.out.print("},");
                }else{
                    System.out.print("}");
                }

        }

        System.out.println("]}");


        
    }

public static void main(String[] args) throws IOException, Arguments.ArgumentException {
    Arguments arguments = new Arguments(
        new Arguments.Option[]{
            new Arguments.RealOption("rate", "Rate parameter"),
            new Arguments.RealOption("bias", "Bias parameter"),
            new Arguments.RealOption("length", "branchlength to evaluate"),
            new Arguments.Option("help", "option to print this message"),
        }
    );


    try {
        arguments.parseArguments(args);
    } catch (Arguments.ArgumentException ae) {
        System.err.println(ae);
        System.exit(1);
    }

    double rate = arguments.getRealOption("rate");
    double bias = arguments.getRealOption("bias");
    double length = arguments.getRealOption("length");

    MarkovRewardDensityCalculator calculator = new MarkovRewardDensityCalculator(rate, bias);

    calculator.run(length);

   
}


private class Result{
    double pdf;
    double cdf;
    double length;
    double s;
    double rate;
    double bias;
    String calculator;


public Result(double pdf, double cdf, double length, double s, double rate, double bias, String calculator) {
    this.pdf = pdf;
    this.cdf = cdf;
    this.length = length;
    this.s = s;
    this.rate = rate;
    this.bias = bias;
    this.calculator = calculator;
}
}
    
}
