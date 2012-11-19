package dr.evomodel.epidemiology.casetocase;

import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.inference.model.AbstractModel;

/**
 * Created with IntelliJ IDEA.
 * User: Matthew Hall
 * Date: 11/05/2012
 * Time: 14:28
 * To change this template use File | Settings | File Templates.
 */


public abstract class AbstractCase extends AbstractModel {

    public AbstractCase(String name){
        super(name);
    }

    protected String caseID;
    //These can be either dates or probability distributions
    protected Object infectionDate;
    protected Object infectiousDate;
    protected Object endOfInfectiousDate;

    public String getName(){
        return caseID;
    }

    public abstract Date getLatestPossibleInfectionDate();

    public abstract Taxa getAssociatedTaxa();

    public abstract boolean culledYet(int time);

    public String toString(){
        return caseID;
    }


}

