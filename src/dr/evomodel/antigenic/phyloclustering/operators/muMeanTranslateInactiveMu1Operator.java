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

public class muMeanTranslateInactiveMu1Operator extends AbstractCoercableOperator   {

	
    private MatrixParameter mu = null;
    private Parameter muMean = null;   

    private Parameter indicators;
    
	
    private double windowSize = 0.5;

	
	public muMeanTranslateInactiveMu1Operator(double weight,  MatrixParameter mu, Parameter indicators, Parameter muMean, double windowSize){
  
        super(CoercionMode.COERCION_ON);

		setWeight(weight);
        this.mu = mu;
        this.indicators = indicators;
        this.muMean = muMean;
		this.windowSize = windowSize;


	}
	
	

	public double doOperation() {
       System.out.println("run here stop");
       System.exit(0);
        //unbounded walk
        double change = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
        
		//change mu1Scale

		double original_muMean_Val = muMean.getParameterValue(0);
   	    double new_muMean_Val = change + original_muMean_Val;
		muMean.setParameterValue(0, new_muMean_Val);
		
		//translate all the inactive mean mu
		int numNodes = mu.getColumnDimension();
		//make sure all the active mu's first dimension stays intact 
		for(int i=0; i < numNodes; i++){
			if( (int) indicators.getParameterValue(i) == 0){
				double oldValue = mu.getParameter(i).getParameterValue(0);
				double newValue =  oldValue +change;
				mu.getParameter(i).setParameterValue(0, newValue);
			}
		}
		
		return 0;

	}
	
	
	


	 //MCMCOperator INTERFACE
   public double getCoercableParameter() {
       return Math.log(windowSize);
   }

   public void setCoercableParameter(double value) {
       windowSize = Math.exp(value);
   }

   public double getRawParameter() {
       return windowSize;
   }

   public double getTargetAcceptanceProbability() {
       return 0.234;
   }

   public double getMinimumAcceptanceLevel() {
       return 0.1;
   }

   public double getMaximumAcceptanceLevel() {
       return 0.4;
   }

   public double getMinimumGoodAcceptanceLevel() {
       return 0.20;
   }

   public double getMaximumGoodAcceptanceLevel() {
       return 0.30;
   }

   public final String getPerformanceSuggestion() {

       double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
       double targetProb = getTargetAcceptanceProbability();

       double ws = OperatorUtils.optimizeWindowSize(windowSize, prob, targetProb);

       if (prob < getMinimumGoodAcceptanceLevel()) {
           return "Try decreasing windowSize to about " + ws;
       } else if (prob > getMaximumGoodAcceptanceLevel()) {
           return "Try increasing windowSize to about " + ws;
       } else return "";
   }

	
	
    public final static String muMeanTranslateInactiveMu1OperatorStr = "muMeanTranslateInactiveMu1Operator";

    public final String getOperatorName() {
        return muMeanTranslateInactiveMu1OperatorStr;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	


    	public final static String  MU = "mu";
    	public final static String  MUMEAN = "muMean";       
    	public final static String INDICATORS = "indicators";
    	public final static String WINDOWSIZE = "windowSize";

        public String getParserName() {
            return muMeanTranslateInactiveMu1OperatorStr;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOWSIZE);
            
                XMLObject cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                cxo = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(MUMEAN);
                Parameter muMean = (Parameter) cxo.getChild(Parameter.class);

            return new muMeanTranslateInactiveMu1Operator(weight, mu, indicators, muMean, windowSize);
        	


        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "changes mu1Scale and make sure the first dimension of the active drifted mus stay the same";
        }

        public Class getReturnType() {
            return muMeanTranslateInactiveMu1Operator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(WINDOWSIZE),
                new ElementRule(MU, Parameter.class),
               new ElementRule(INDICATORS, Parameter.class),
               new ElementRule(MUMEAN, Parameter.class),

        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
