package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class SmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable {

    public SmoothSkygridLikelihood(String name, IntervalList intervalList) {
        super(name, intervalList);
    }

    @Override
    public Type getUnits() {
        return null;
    }

    @Override
    public void setUnits(Type units) {

    }

    @Override
    protected double calculateLogLikelihood() {
        return 0;
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        return 0;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        return 0;
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }
}
