package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.NodeRef;
import dr.inference.model.AbstractModel;
import dr.inference.model.Variable;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: Matthew Hall
 * Date: 22/05/2012
 * Time: 14:40
 * To change this template use File | Settings | File Templates.
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

    public abstract double P(AbstractCase[] paintings, double[] times);

    public abstract double logP(AbstractCase[] paintings, double[] times);

    public abstract double probInfectiousBy(AbstractCase painting, double time);

    public abstract double logProbInfectiousBy(AbstractCase painting, double time);

}
