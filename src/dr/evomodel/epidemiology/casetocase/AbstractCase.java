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
    protected Date examDate;
    protected Date endOfInfectiousDate;

    public String getName(){
        return caseID;
    }

    public Taxa getAssociatedTaxa() {
        return associatedTaxa;
    }

    public Date getExamDate() {
        return examDate;
    }

    public double getExamTime() {
        return examDate.getTimeValue();
    }

    public Date getCullDate(){
        return endOfInfectiousDate;
    }

    public double getCullTime(){
        return endOfInfectiousDate.getTimeValue();
    }

    public abstract boolean culledYet(double time);

    public String toString(){
        return caseID;
    }

    public abstract double[] getCoords();


}

