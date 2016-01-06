package dr.inference.operators;

import dr.math.distributions.NormalDistribution;

/**
 * Created by max on 12/3/15.
 */
public abstract class AbstractHamiltonianMCOperator extends AbstractCoercableOperator{
    public AbstractHamiltonianMCOperator(CoercionMode mode, double momentumSd) {
        super(mode);
        this.momentumSd=momentumSd;
    }

    protected double getMomentumSd()
    {return momentumSd;}

    protected void setMomentumSd(double momentum){
        momentumSd=momentum;
    }

    private double momentumSd;
    protected double[] momentum;

    protected void drawMomentum(int size){
        momentum=new double[size];
        for (int i = 0; i <size ; i++) {
            momentum[i]= (Double) (new NormalDistribution(0.0, momentumSd)).nextRandom();
        }
    }
}
