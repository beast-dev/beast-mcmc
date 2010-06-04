package dr.evomodel.graph;

import java.util.ArrayList;
import java.util.List;

import dr.inference.model.CompoundParameter;
import dr.inference.model.IntersectionBounds;
import dr.inference.model.Parameter;

/*
 * Implements a version of the CompoundParameter that can do store/restore 
 * after adding and removing variables
 */
public class VariableDimensionCompoundParameter extends CompoundParameter {

	public VariableDimensionCompoundParameter(Parameter[] params) {
		super(params);
	}

	public VariableDimensionCompoundParameter(String name) {
		super(name);
	}

    protected void storeValues() {
    	super.storeValues();
    	copyValues(parameters, storedParameters);
    	copyValues(uniqueParameters, storedUniqueParameters);
    	copyValues(pindex, storedPindex);
    	if(bounds!=null){
//			storedBounds = new dr.inference.model.IntersectionBounds();
	    	storedBounds = bounds;
	    	System.err.println("fixme!");
    	}
    	storedDimension = dimension;
    }
    
    protected void restoreValues() {
    	copyValues(storedParameters, parameters);
    	copyValues(storedUniqueParameters, uniqueParameters);
    	copyValues(storedPindex, pindex);
    	dimension = storedDimension;
    	super.restoreValues();
    }
    protected void acceptValues() {
    }

    private <E> void copyValues(List<E> from, List<E> to){
    	to.clear();
    	to.addAll(from);
    }

    protected final List<Parameter> storedUniqueParameters = new ArrayList<Parameter>();
    protected final ArrayList<Parameter> storedParameters = new ArrayList<Parameter>();
    protected final ArrayList<Integer> storedPindex = new ArrayList<Integer>();
    protected IntersectionBounds storedBounds = null;
    protected int storedDimension;

}
