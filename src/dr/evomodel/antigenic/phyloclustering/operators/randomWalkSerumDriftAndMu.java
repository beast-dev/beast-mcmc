package dr.evomodel.antigenic.phyloclustering.operators;


import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;


/**
 * An operator to cluster viruses using a phylogenetic tree
 *
 * @author Charles Cheung
 * @author Trevor Bedford
 */
public class randomWalkSerumDriftAndMu extends SimpleMCMCOperator  {
	
	
    public static final String SERUMDRIFT_AND_MU_OPERATOR = "serumDriftAndMuOperator";
    
    
	//Variables
    Parameter indicators;
    MatrixParameter mu; //mu - means
    Parameter serumDrift;	
    private TreeModel treeModel;

    private double maxWalkSize;
	 
    //Constructor
    public randomWalkSerumDriftAndMu(	MatrixParameter mu, 
    									double weight, 
    									Parameter indicatorsParameter, 
    									Parameter serumDrift_in,
    									double max_walk_size_in,
    									TreeModel treeModel_in) {
    	this.mu = mu;
    	this.indicators = indicatorsParameter;
       	this.serumDrift = serumDrift_in;
       	this.maxWalkSize = max_walk_size_in;
		this.treeModel= treeModel_in;

        setWeight(weight);
        
        System.out.println("Finished loading the constructor for SERUMDRIFT_AND_MU_OPERATOR");
    	
    }
    

    
 
    /**
     * change the parameter and return the log hastings ratio.
     */
    public final double doOperation() {

    	double logHastingRatio = 0; //initiate the log Metropolis Hastings ratio of the MCMC

    	
    	int rootNode  = treeModel.getRoot().getNumber();
        //perform proposal
    	
    	//random walk serum drift 1
		double change = MathUtils.nextDouble()*maxWalkSize- maxWalkSize/2 ;
		double originalValue = serumDrift.getParameterValue(0);
		double newValue = originalValue + change;
		serumDrift.setParameterValue(0, newValue);
		
		for(int i=0; i < mu.getParameterCount(); i++){
			if( (int) indicators.getParameterValue(i) == 1  && i != rootNode ){
				Parameter curMu = mu.getParameter(i);
				double originalMu0 = curMu.getParameterValue(0);
				double newMu0 = originalMu0 * newValue/originalValue;
				curMu.setParameterValue(0, newMu0);
			}
		}
    	
			
    	return(logHastingRatio);    	
    }
    	

	  
	public void accept(double deviation) {
    	super.accept(deviation);         	
    }
    
    public void reject(){
    	super.reject();
    }
    
	

             
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return SERUMDRIFT_AND_MU_OPERATOR;
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

            public String getPerformanceSuggestion() {
                if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
                    return "";
                } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
                    return "";
                } else {
                    return "";
                }
            }

        
           
        

            public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
            	

            	public final static String  MU = "mu";
            	public final static String SERUMDRIFT = "serumDrift";
            	public final static String INDICATORS = "indicators";
            	public final static String WALKSIZE = "walkSize";
                
                public String getParserName() {
                    return SERUMDRIFT_AND_MU_OPERATOR;
                }

                /* (non-Javadoc)
                 * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
                 */
                public Object parseXMLObject(XMLObject xo) throws XMLParseException {

                    double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

                    double walk_size = 0.05;
                	if (xo.hasAttribute(WALKSIZE)) {
                		walk_size = xo.getDoubleAttribute(WALKSIZE);
                	}
                    
                    XMLObject cxo =  xo.getChild(MU);
                        MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        cxo = xo.getChild(SERUMDRIFT);
                        Parameter serumDrift = (Parameter) cxo.getChild(Parameter.class);

                        cxo = xo.getChild(INDICATORS);
                        Parameter indicators = (Parameter) cxo.getChild(Parameter.class);

                        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);


                        return new randomWalkSerumDriftAndMu(mu,  weight,  indicators, serumDrift, walk_size, treeModel);
                    

                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
                }

                public Class getReturnType() {
                    return TreeClusterAlgorithmOperator.class;
                }


                public XMLSyntaxRule[] getSyntaxRules() {
                    return rules;
                }

                private final XMLSyntaxRule[] rules = {
                        AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                        AttributeRule.newDoubleRule(WALKSIZE),
                        new ElementRule(MU, Parameter.class),
                        new ElementRule(SERUMDRIFT, Parameter.class),
                       new ElementRule(INDICATORS, Parameter.class),
                       new ElementRule(TreeModel.class),
                };
            
            };


        }



        
