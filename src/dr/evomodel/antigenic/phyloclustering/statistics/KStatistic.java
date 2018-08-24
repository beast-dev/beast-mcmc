



package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.inference.model.*;
import dr.xml.*;

/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */



public class KStatistic extends Statistic.Abstract implements VariableListener {

	
		
    public static final String K_STATISTIC = "kStatistic";

    public KStatistic(Parameter indicators) {
        this.indicatorsParameter = indicators;
        indicatorsParameter.addParameterListener(this);
    }
    


    public int getDimension() {
        return 1;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {

    	double count = 0;
   		for(int i=0; i < indicatorsParameter.getDimension(); i++){
   			if( (int) indicatorsParameter.getParameterValue(i) == 1 ){
   				count++;
   			}
   		}

       return count;

    }

    public String getDimensionName(int dim) {
    	String name = "K";
        return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String INDICATORS = "indicators";

        public String getParserName() {
            return K_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS);
            return new KStatistic(indicators);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return KStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(INDICATORS, Parameter.class)
        };
    };

    private Parameter indicatorsParameter;

}
