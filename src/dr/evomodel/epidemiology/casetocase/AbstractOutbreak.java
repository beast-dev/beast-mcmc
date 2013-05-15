package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.AbstractModel;

import java.util.ArrayList;

/**
 * Abstract class for outbreaks
 *
 * User: Matthew Hall
 * Date: 14/04/13
 */

public abstract class AbstractOutbreak extends AbstractModel {

    public AbstractOutbreak(String name){
        super(name);
    }

    protected ArrayList<AbstractCase> cases;

    public ArrayList<AbstractCase> getCases(){
        return new ArrayList<AbstractCase>(cases);
    }

    public int size(){
        return cases.size();
    }

    public AbstractCase getCase(int i){
        return cases.get(i);
    }

    public AbstractCase getCase(String name){
        for(AbstractCase thisCase: cases){
            if(thisCase.getName().equals(name)){
                return thisCase;
            }
        }
        return null;
    }

    public abstract double P(AbstractCase parentPainting, AbstractCase childPainting, double time1,
                             double time2, boolean extended);

    public abstract double logP(AbstractCase parentPainting, AbstractCase childPainting, double time1,
                             double time2, boolean extended);

    public abstract double probInfectiousBy(AbstractCase painting, double time, boolean extended);

    public abstract double logProbInfectiousBy(AbstractCase painting, double time, boolean extended);

}
