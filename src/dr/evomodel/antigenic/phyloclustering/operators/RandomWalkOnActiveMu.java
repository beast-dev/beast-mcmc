
package dr.evomodel.antigenic.phyloclustering.operators;

import java.util.LinkedList;

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

public class RandomWalkOnActiveMu extends AbstractCoercableOperator {

	
    private MatrixParameter mu = null;
    private MatrixParameter virusLocations = null;
    private MatrixParameter virusLocationsTreeNode = null;
    private Parameter indicators;
    private TreeModel treeModel;
    
    private int numdata;   
    private int numNodes;
	private int []correspondingTreeIndexForVirus = null; //relates treeModels's indexing system to cluster label's indexing system of viruses. Gets assigned

    private double windowSize = 0.5;

	
	public RandomWalkOnActiveMu(double weight, MatrixParameter virusLocations, MatrixParameter mu, Parameter indicators,  TreeModel treeModel_in, double windowSize, MatrixParameter virusLocationsTreeNode_in){
    
        super(CoercionMode.COERCION_ON);
		
		setWeight(weight);
        this.windowSize = windowSize;

        this.virusLocations = virusLocations;
        this.mu = mu;
        this.indicators = indicators;
		this.treeModel= treeModel_in;
		this.virusLocationsTreeNode = virusLocationsTreeNode_in;
		
		numNodes = treeModel.getNodeCount();
		numdata = virusLocations.getColumnDimension();

		correspondingTreeIndexForVirus = TreeClusteringSharedRoutines.setMembershipTreeToVirusIndexes(numdata, virusLocations, numNodes, treeModel);
    	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
	}
	
	

	public double doOperation() {


        //first, randomly select an "on" node to overwrite
		int originalNode = TreeClusteringSharedRoutines.findAnOnNodeIncludingRootRandomly(numNodes, indicators);			//find an on-node
		//unbounded walk
		int dimSelect = (int) Math.floor(  MathUtils.nextDouble()* 2 );
        double change = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
		double originalValue = mu.getParameter(originalNode).getParameterValue(dimSelect);
		mu.getParameter(originalNode ).setParameterValue(dimSelect, originalValue + change);
	
		//a. by removing the selected node, each child of this node should be updated to keep the absolute location of 
		//the child cluster fixed as before
		LinkedList<Integer> childrenOriginalNode = TreeClusteringSharedRoutines.findActiveBreakpointsChildren(originalNode, numNodes, treeModel, indicators);
		for(int i=0; i < childrenOriginalNode.size(); i++){
			int muIndexNum = childrenOriginalNode.get(i).intValue() ;
			Parameter curMu = mu.getParameter( muIndexNum );
			double curMu_original = curMu.getParameterValue( dimSelect);
			mu.getParameter(muIndexNum).setParameterValue(dimSelect, curMu_original - change);
		}
		
		
		
		
		
		//the virus location needs to be updated because the mu's are updated 	  				
    	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
		

        return 0.0;
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

    
    
    
    public final static String RANDOMWALKACTIVEMU = "randomWalkOnActiveMu";

    public final String getOperatorName() {
        return RANDOMWALKACTIVEMU;
    }

    


    
    
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	

        public final static String VIRUSLOCATIONS = "virusLocations";
        public final static String VIRUSLOCATIONSTREENODE = "virusLocationsTreeNodes";
    	public final static String  MU = "mu";
    	public final static String INDICATORS = "indicators";
    	public final static String WINDOWSIZE = "windowSize";


        public String getParserName() {
            return RANDOMWALKACTIVEMU;
        }

        /* (non-Javadoc)
         * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOWSIZE);

            
            XMLObject cxo = xo.getChild(VIRUSLOCATIONS);
                MatrixParameter virusLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
                cxo = xo.getChild(VIRUSLOCATIONSTREENODE);
                MatrixParameter virusLocationsTreeNode = (MatrixParameter) cxo.getChild(MatrixParameter.class);
               
                
                cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                cxo = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

                TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            return new RandomWalkOnActiveMu(weight, virusLocations, mu, indicators, treeModel, windowSize, virusLocationsTreeNode);
            

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "changes serum drift and make sure the first dimension of the active drifted mus stay the same";
        }

        public Class getReturnType() {
            return RandomWalkOnActiveMu.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(WINDOWSIZE),
                new ElementRule(VIRUSLOCATIONS, Parameter.class),
                new ElementRule(VIRUSLOCATIONSTREENODE, MatrixParameter.class),
                new ElementRule(MU, Parameter.class),
               new ElementRule(INDICATORS, Parameter.class),
               new ElementRule(TreeModel.class),

        };
    
    };



    public int getStepCount() {
        return 1;
    }
    

}
