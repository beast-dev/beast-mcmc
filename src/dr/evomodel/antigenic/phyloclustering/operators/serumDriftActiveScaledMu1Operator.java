
package dr.evomodel.antigenic.phyloclustering.operators;

import dr.evomodel.antigenic.phyloclustering.TreeClusteringSharedRoutines;
import dr.evomodel.tree.TreeModel;
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

public class serumDriftActiveScaledMu1Operator extends AbstractCoercableOperator {

	
    private MatrixParameter mu = null;
    private Parameter serumDrift = null;   
    private MatrixParameter virusLocations = null;
    private MatrixParameter virusLocationsTreeNode = null;
    private Parameter indicators;
    private TreeModel treeModel;
    
    private int numdata;   
    private int numNodes;
	private int []correspondingTreeIndexForVirus = null; //relates treeModels's indexing system to cluster label's indexing system of viruses. Gets assigned

	private double scaleFactor;

	
	public serumDriftActiveScaledMu1Operator(double weight, MatrixParameter virusLocations, MatrixParameter mu, Parameter indicators, Parameter serumDrift, TreeModel treeModel_in, double scale, MatrixParameter virusLocationsTreeNode_in){
    
        super(CoercionMode.COERCION_ON);
		
		setWeight(weight);
        this.virusLocations = virusLocations;
        this.mu = mu;
        this.indicators = indicators;
        this.serumDrift = serumDrift;
		this.treeModel= treeModel_in;
		this.scaleFactor = scale;
		this.virusLocationsTreeNode = virusLocationsTreeNode_in;
		
		numNodes = treeModel.getNodeCount();
		numdata = virusLocations.getColumnDimension();

		correspondingTreeIndexForVirus = TreeClusteringSharedRoutines.setMembershipTreeToVirusIndexes(numdata, virusLocations, numNodes, treeModel);
    	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
	}
	
	

	public double doOperation() {

		
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

 
		//change serum drift
        //double WALK_SIZE = 0.2; //when the walk size becomes 0.5, .... things become weird.. too big
		//double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ;
		double original_serumDrift_Val = serumDrift.getParameterValue(0);
	//	System.out.println("original_serumDrift_Val = " + original_serumDrift_Val);
		//double new_serumDrift_Val = change + original_serumDrift_Val;
	   double new_serumDrift_Val = scale * original_serumDrift_Val;
	   
	   
		serumDrift.setParameterValue(0, new_serumDrift_Val);
		
	//	System.out.println("new_serumDrift_Val=" + serumDrift.getParameterValue(0));
		//make sure all the active mu's first dimension stays intact 
		for(int i=0; i < numNodes; i++){
			if( (int) indicators.getParameterValue(i) == 1){
				double oldValue = mu.getParameter(i).getParameterValue(0);				
				double newValue =  oldValue * new_serumDrift_Val/original_serumDrift_Val;
				mu.getParameter(i).setParameterValue(0, newValue);
	//			System.out.println("indicator" + i + "  oldValue = " + oldValue + " and newValue=" + mu.getParameter(i).getParameterValue(0));
			}
		}
		
		//the virus location needs to be updated because the mu's are updated 	  				
    	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
		
		
        double logq = -Math.log(scale);
    //    System.out.println("logq=" + logq);
    //    System.out.println("================================================");
		//return 0;
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
	
	
/*	
    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
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
*/
    
    
    
    public final static String SERUMDRIFTACTIVESCALEDMU1Operator = "serumDriftActiveScaledMu1Operator";

    public final String getOperatorName() {
        return SERUMDRIFTACTIVESCALEDMU1Operator;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	

        public final static String VIRUSLOCATIONS = "virusLocations";
        public final static String VIRUSLOCATIONSTREENODE = "virusLocationsTreeNodes";
    	public final static String  MU = "mu";
    	public final static String  SERUMDRIFT = "serumDrift";       
    	public final static String INDICATORS = "indicators";
    	public final static String SCALE = "scaleFactor";


        public String getParserName() {
            return SERUMDRIFTACTIVESCALEDMU1Operator;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double scale = xo.getDoubleAttribute(SCALE);

            
            XMLObject cxo = xo.getChild(VIRUSLOCATIONS);
                MatrixParameter virusLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
                cxo = xo.getChild(VIRUSLOCATIONSTREENODE);
                MatrixParameter virusLocationsTreeNode = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
                
                cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                cxo = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(SERUMDRIFT);
                Parameter serumDrift = (Parameter) cxo.getChild(Parameter.class);

                TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            return new serumDriftActiveScaledMu1Operator(weight, virusLocations, mu, indicators, serumDrift, treeModel, scale, virusLocationsTreeNode);
            

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
                new ElementRule(VIRUSLOCATIONS, Parameter.class),
                new ElementRule(VIRUSLOCATIONSTREENODE, MatrixParameter.class),
                new ElementRule(MU, Parameter.class),
               new ElementRule(INDICATORS, Parameter.class),
               new ElementRule(SERUMDRIFT, Parameter.class),
               new ElementRule(TreeModel.class),

        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
