
package dr.evomodel.antigenic.phyloclustering.statistics;

import java.util.LinkedList;

import dr.inference.model.*;
import dr.xml.*;

/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */

public class ActiveIndicatorsStatistic extends Statistic.Abstract implements VariableListener {

	private LinkedList<Double> activeNodes  = new LinkedList<Double>();
	
	private int max_dim;
	
    public static final String ACTIVE_INDICATORS_STATISTIC = "activeIndicatorsStatistic";

    public ActiveIndicatorsStatistic(Parameter indicators, int maxDim_in) {
        this.indicatorsParameter = indicators;
        indicatorsParameter.addParameterListener(this);
        max_dim = maxDim_in;
    }
    


    public int getDimension() {
        return max_dim;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {

    	if(dim ==0){
    		activeNodes  = new LinkedList<Double>();  // reset linkedlist	
    	    //determine all the nodes that are active.
    		for(int i=0; i < indicatorsParameter.getDimension(); i++){
    			if( (int) indicatorsParameter.getParameterValue(i) == 1 ){
    				activeNodes.addLast(new Double(i));
    			}
    		}
    		//System.out.println("active node size is = " + activeNodes.size() );
    	}
    	
    	
        double val = -1;
        if(dim < activeNodes.size()){
        	val = activeNodes.get(dim).doubleValue();
        }
        
        
        //if the number of active nodes is more than the max number allowed, then for the last placeholder, print -9999999 instead
        //to indicate that there are more than MAX_DIM of active nodes and we aren't able to print them all.
        if(dim == (max_dim -1) && ( activeNodes.size() > max_dim) ){
        	val = -9999999;
        }

       return val;

    }

    public String getDimensionName(int dim) {
    	String name = "on_" + (dim+1);
        return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String INDICATORS = "indicators";
        public final static String MAXDIMSTR = "maxDim";


        public String getParserName() {
            return ACTIVE_INDICATORS_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS);
            
        	int maxDim = 30;
        	if(xo.hasAttribute(MAXDIMSTR)){
        		maxDim = xo.getIntegerAttribute(MAXDIMSTR);
        	}

            return new ActiveIndicatorsStatistic(indicators, maxDim);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return ActiveIndicatorsStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(INDICATORS, Parameter.class),
            AttributeRule.newDoubleRule(MAXDIMSTR, true, "the variance of mu"),

        };
    };

    private Parameter indicatorsParameter;

}
