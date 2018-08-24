
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

public class muPrecisionInactiveMuOperator extends AbstractCoercableOperator {

	
   
    private Parameter muPrecision = null;   
    private MatrixParameter mu = null;
    private Parameter indicators = null;
    private double scaleFactor;
    private Parameter muMean = null;

	public muPrecisionInactiveMuOperator(double weight, MatrixParameter mu, Parameter muPrec,  double scale, Parameter indicators, Parameter muMean){
    
        super(CoercionMode.COERCION_ON);
		
		setWeight(weight);
        this.mu = mu;
        this.muPrecision = muPrec;
		this.scaleFactor = scale;
		this.indicators = indicators;
		this.muMean = muMean;
	
	}
	
	

	public double doOperation() {

		
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
 
		double original_muPrec_Val = muPrecision.getParameterValue(0);
        double new_muPrec_Val = scale * original_muPrec_Val;
	   
	   
		muPrecision.setParameterValue(0, new_muPrec_Val);
		 
		
		for(int i=0; i < mu.getColumnDimension(); i++){
			
			//only change the inactive mus
			if(indicators.getParameterValue(i) == 0){
				for(int j=0; j < 2; j++){
					if(j==0){
						double translatedValue = (mu.getParameter(i).getParameterValue(j) - muMean.getParameterValue(0) );				
						double newValue =  translatedValue *  Math.sqrt( original_muPrec_Val/new_muPrec_Val) + muMean.getParameterValue(0); 
						mu.getParameter(i).setParameterValue(j, newValue);
					}
					else{
						double oldValue = mu.getParameter(i).getParameterValue(j);				
						double newValue =  oldValue *  Math.sqrt( original_muPrec_Val/new_muPrec_Val); 
						mu.getParameter(i).setParameterValue(j, newValue);
					}
				}
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
	
	

    
    public final static String CLASSNAME = "muPrecisionInactiveMuOperator";

    public final String getOperatorName() {
        return CLASSNAME;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	

        public final static String MU = "mu";
    	public final static String  MUPREC = "muPrec";       

    	public final static String SCALE = "scaleFactor";

    	public final static String MUMEAN = "muMean";
    	public final static String INDICATORS = "indicators";


        public String getParserName() {
            return CLASSNAME;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double scale = xo.getDoubleAttribute(SCALE);

            
            XMLObject cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
    
                cxo = xo.getChild(MUPREC);
                Parameter muPrec = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
                
                cxo = xo.getChild(MUMEAN);
                Parameter muMean = (Parameter) cxo.getChild(Parameter.class);
                
            return new muPrecisionInactiveMuOperator(weight, mu, muPrec,  scale, indicators, muMean);
            

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "changes serum drift and make sure the first dimension of the active drifted mus stay the same";
        }

        public Class getReturnType() {
            return muPrecisionInactiveMuOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(MU, Parameter.class),
               new ElementRule(MUPREC, Parameter.class),
               new ElementRule(INDICATORS, Parameter.class),
               new ElementRule(MUMEAN, Parameter.class)
        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
