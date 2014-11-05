package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.inference.model.AbstractModel;

/**
 * Abstract class for outbreak; best implemented as an inner class in implementations of AbstractOutbreak
 *
 * User: Matthew Hall
 * Date: 15/04/13
 */


public abstract class AbstractCase extends AbstractModel {

    public AbstractCase(String name){
        super(name);
    }

    protected String caseID;
    protected Taxa associatedTaxa;
    protected double examTime;
    protected double endOfInfectiousTime;
    protected boolean wasEverInfected;

    public String getName(){
        return caseID;
    }

    public Taxa getAssociatedTaxa() {
        return associatedTaxa;
    }

    public double getExamTime() {
        return examTime;
    }

    public double getCullTime(){
        return endOfInfectiousTime;
    }

    public abstract boolean culledYet(double time);

    public String toString(){
        return caseID;
    }

    public abstract double[] getCoords();

    public boolean wasEverInfected(){
        return wasEverInfected;
    }

    public void setEverInfected(boolean value){
        wasEverInfected = value;
    }


}

