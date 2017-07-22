package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.inference.model.Parameter;
import dr.inference.model.*;
import dr.xml.*;


/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */

public class indicatorsStatistic  extends Statistic.Abstract implements VariableListener {
	
    private Parameter indicators;
    
    public static final String INDICATORS_STATISTIC = "indicatorsStatistic";

    public indicatorsStatistic( Parameter indicators) {
        this.indicators = indicators;
        indicators.addParameterListener(this);
    }
    


    public int getDimension() {
        return indicators.getDimension();
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {
        return (  indicators.getParameterValue(dim)  );

    }

    
    
     
    
    
    public String getDimensionName(int dim) {
    	String name = "indicators_" +  dim ;
        return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String INDICATORS_STRING = "indicators";

        public String getParserName() {
            return INDICATORS_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS_STRING);

            return new indicatorsStatistic( indicators);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return indicatorsStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(INDICATORS_STRING, Parameter.class)

        };
    };

    

}
