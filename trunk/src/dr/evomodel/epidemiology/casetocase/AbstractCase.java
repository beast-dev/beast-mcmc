package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.inference.model.AbstractModel;

/**
 * Abstract class for cases; best implemented as an inner class in implementations of AbstractOutbreak
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


    public String getName(){
        return caseID;
    }

    public abstract Date getLatestPossibleInfectionDate();

    public Taxa getAssociatedTaxa() {
        return associatedTaxa;
    }

    public abstract boolean culledYet(int time);

    public String toString(){
        return caseID;
    }




}

