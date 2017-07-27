
package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import dr.evomodel.antigenic.phyloclustering.operators.serumDriftActiveScaledMu1Operator;
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

public class serumDriftScalingMuOperator extends AbstractCoercableOperator {

	

    private Parameter serumDrift = null;   

    private MatrixParameter mu = null;
    private Parameter muMean = null;
    private Parameter muPrec = null;
	private double scaleFactor;

	
	public serumDriftScalingMuOperator(double weight, MatrixParameter mu, Parameter muMean, Parameter muPrec, Parameter serumDrift, double scale){
    
        super(CoercionMode.COERCION_ON);
		
		setWeight(weight);
        this.mu = mu;
        this.muMean = muMean;
        this.muPrec = muPrec;
        this.serumDrift = serumDrift;
		this.scaleFactor = scale;
	}
	
	

	public double doOperation() {

		
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

     //   System.out.println("serumDriftScaling operator ran");
      //  System.out.println("scale=" + scale);
        //changing serum drift
        serumDrift.setParameterValue(0, scale *serumDrift.getParameterValue(0) );
		
		//changing mu
		//System.out.println("dimension=" + mu.getColumnDimension());
		for(int i=0; i < mu.getColumnDimension(); i++){
			Parameter m = mu.getParameter(i);
			m.setParameterValue(0, scale*m.getParameterValue(0));
		}
		//changing muMean
		muMean.setParameterValue(0, scale*muMean.getParameterValue(0));
		
		//changing muPrec
		muPrec.setParameterValue(0, scale*scale*muPrec.getParameterValue(0));
				
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
	
	

    
    public final static String SERUMDRIFTSCALINGMUOperator = "serumDriftScalingMuOperator";

    public final String getOperatorName() {
        return SERUMDRIFTSCALINGMUOperator;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	

    	public final static String SCALE = "scaleFactor";
    	public final static String  MU = "mu";
    	public final static String  SERUMDRIFT = "serumDrift";       
    	public final static String MUMEAN = "muMean";
    	public final static String MUPREC = "muPrec";

        public String getParserName() {
            return SERUMDRIFTSCALINGMUOperator;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double scale = xo.getDoubleAttribute(SCALE);              
                
            XMLObject cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                cxo = xo.getChild(SERUMDRIFT);
                Parameter serumDrift = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(MUMEAN);
                Parameter muMean = (Parameter) cxo.getChild(Parameter.class);
                
                
                cxo = xo.getChild(MUPREC);
                Parameter muPrec = (Parameter) cxo.getChild(Parameter.class);
                
                
            return new serumDriftScalingMuOperator(weight, mu, muMean, muPrec, serumDrift,  scale);
            
            //	public serumDriftScalingMuOperator(double weight, MatrixParameter mu, Parameter muMean, Parameter muPrec, Parameter serumDrift, double scale){


        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "changes serum drift and make sure the first dimension of the active drifted mus stay the same";
        }

        public Class getReturnType() {
            return serumDriftActiveScaledMu1Operator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(SCALE),
                new ElementRule(MU, Parameter.class),
                new ElementRule(MUMEAN, Parameter.class),
                new ElementRule(MUPREC, Parameter.class),
                new ElementRule(SERUMDRIFT, Parameter.class),
        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
