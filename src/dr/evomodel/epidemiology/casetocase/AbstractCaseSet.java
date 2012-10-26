package dr.evomodel.epidemiology.casetocase;

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
public abstract class AbstractCaseSet extends AbstractModel {

    public AbstractCaseSet(String name){
        super(name);
    }

    protected ArrayList<AbstractCase> cases;

    /*Likelihood for the root if it is painted with 'farm'  */

    public abstract double rootBranchLikelihood(AbstractCase farm, Integer farmInfectiousBy);

    public abstract double rootBranchLogLikelihood(AbstractCase farm, Integer farmInfectiousBy);


    /*Likelihood for a branch leading from a node painted 'parent' to a node painted 'child'. This deals with the
    probability that farm2 is infected at the first time and infectious at the second. The calculations for farm1
    will be done when that node is considered as a child, or as the root*/

    public abstract double branchLikelihood(AbstractCase parent, AbstractCase child, Integer childInfected,
                                            Integer childInfectiousBy);

    public abstract double branchLogLikelihood(AbstractCase parent, AbstractCase child, Integer childInfected,
                                            Integer childInfectiousBy);

    public abstract ArrayList<AbstractCase> getCases();

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

}
