
package dr.evomodel.antigenic.phyloclustering.operators;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class serumPrecisionSerumLocOperator extends AbstractCoercableOperator {

	
   
    private Parameter serumPrecision = null;   
    private MatrixParameter serumLocations = null;

    private double scaleFactor;

	public serumPrecisionSerumLocOperator(double weight, MatrixParameter serumLocations, Parameter serumPrec,  double scale){
    
        super(CoercionMode.COERCION_ON);
		
		setWeight(weight);
        this.serumLocations = serumLocations;
        this.serumPrecision = serumPrec;
		this.scaleFactor = scale;
	
	}
	
	

	public double doOperation() {

		
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
 
		double original_serumPrec_Val = serumPrecision.getParameterValue(0);
        double new_serumPrec_Val = scale * original_serumPrec_Val;
	   
	   
		serumPrecision.setParameterValue(0, new_serumPrec_Val);
		 
		
		for(int i=0; i < serumLocations.getColumnDimension(); i++){
				for(int j=0; j < 2; j++){
					double oldValue = serumLocations.getParameter(i).getParameterValue(j);				
					double newValue =  oldValue *  Math.sqrt( original_serumPrec_Val/new_serumPrec_Val); 
					serumLocations.getParameter(i).setParameterValue(j, newValue);
				}

		}
		
		
        double logq = -Math.log(scale);
        return logq;
	}
	
	
	


	//copied from the original ScaleOperator
    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

	//copied from the original ScaleOperator
    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

	//copied from the original ScaleOperator
    public double getRawParameter() {
        return scaleFactor;
    }

	
	
	//copied from the original ScaleOperator
    public double getTargetAcceptanceProbability() {
        return 0.234;
    }
	//copied from the original ScaleOperator
    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
	
	

    
    public final static String SERUMPRECSCALEALLSERUMLOC = "serumPrecScaleAllSerumLoc";

    public final String getOperatorName() {
        return SERUMPRECSCALEALLSERUMLOC;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	

        public final static String SERUMLOCATIONS = "serumLocations";
    	public final static String  SERUMPREC = "serumPrec";       

    	public final static String SCALE = "scaleFactor";


        public String getParserName() {
            return SERUMPRECSCALEALLSERUMLOC;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double scale = xo.getDoubleAttribute(SCALE);

            
            XMLObject cxo = xo.getChild(SERUMLOCATIONS);
                MatrixParameter serumLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
    
                cxo = xo.getChild(SERUMPREC);
                Parameter serumPrec = (Parameter) cxo.getChild(Parameter.class);

            return new serumPrecisionSerumLocOperator(weight, serumLocations, serumPrec,  scale);
            

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "changes serum drift and make sure the first dimension of the active drifted mus stay the same";
        }

        public Class getReturnType() {
            return serumPrecisionSerumLocOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(SERUMLOCATIONS, Parameter.class),
               new ElementRule(SERUMPREC, Parameter.class),
        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
