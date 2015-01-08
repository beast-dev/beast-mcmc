package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.prior.Prior;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.util.Attribute;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements a generic univariate slice sampler.
 *
 * See: RM Neal (2003) Slice Sampling, Annals of Statistics, 31, 705-767 (with discussion)
 *
 * @author Marc A. Suchard
 */
public class SliceOperator extends SimpleMetropolizedGibbsOperator {

    public SliceOperator(Variable<Double> variable) {
        this(new SliceInterval.SteppingOut(), variable);
    }

    public SliceOperator(SliceInterval sliceInterval, Variable<Double> variable) {
        this.sliceInterval = sliceInterval;

        if (variable.getSize() != 1) {
            throw new RuntimeException("Generic slice sampler is currently for univariate parameters only");
        }        
        this.variable = variable;
        sliceInterval.setSliceSampler(this);
    }

    public Variable<Double> getVariable() {
        return variable;
    }

    public double doOperation(Prior prior, Likelihood likelihood) throws OperatorFailedException {
        double logPosterior = evaluate(likelihood, prior, 1.0);
        double cutoffDensity = logPosterior + MathUtils.randomLogDouble();
        sliceInterval.drawFromInterval(prior, likelihood, cutoffDensity, width);
        // No need to set variable, as SliceInterval has already done this (and recomputed posterior)
        return 0;
    }

    public int getStepCount() {
        return 1;
    }

    public String getOperatorName() {
        return "genericSliceSampler";
    }

    public static void main(String[] arg) {

        // Define normal model
        Parameter meanParameter = new Parameter.Default(1.0); // Starting value
        Variable<Double> stdev = new Variable.D(1.0, 1); // Fixed value
        ParametricDistributionModel densityModel = new NormalDistributionModel(meanParameter, stdev);
        DistributionLikelihood likelihood = new DistributionLikelihood(densityModel);

        // Define prior
        DistributionLikelihood prior = new DistributionLikelihood(new NormalDistribution(0.0, 1.0)); // Hyper-priors
        prior.addData(meanParameter);

        // Define data
        likelihood.addData(new Attribute.Default<double[]>("Data", new double[] {0.0, 1.0, 2.0}));

        List<Likelihood> list = new ArrayList<Likelihood>();
        list.add(likelihood);
        list.add(prior);
        CompoundLikelihood posterior = new CompoundLikelihood(0, list);
        SliceOperator sliceSampler = new SliceOperator(meanParameter);

        final int length = 10000;
        double mean = 0;
        double variance = 0;

        for(int i = 0; i < length; i++) {
            try {
                sliceSampler.doOperation(null, posterior);
            } catch (OperatorFailedException e) {
                System.err.println(e.getMessage());
            }
            double x = meanParameter.getValue(0);
            mean += x;
            variance += x*x;
        }
        mean /= length;
        variance /= length;
        variance -= mean*mean;
        System.out.println("E(x) = "+mean);
        System.out.println("V(x) = "+variance);
    }

    private final SliceInterval sliceInterval;
    private final double width = 1.0;
    private final Variable<Double> variable;
}
