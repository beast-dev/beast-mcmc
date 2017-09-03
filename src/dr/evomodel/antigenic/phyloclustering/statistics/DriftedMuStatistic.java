
package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */

public class DriftedMuStatistic extends Statistic.Abstract implements VariableListener {

	static  int MAX_DIM = 30;
	

	private MatrixParameter mu;
    private TreeModel treeModel;
    private Parameter indicators;
//    private Parameter locationDrift;
    private Parameter mu1ScaleParameter;
    private Parameter mu2ScaleParameter;
    private Parameter muMeanParameter;
    
    public static final String DRIFTED_MU_STATISTIC = "driftedMuStatistic";

   // public DriftedMuStatistic( TreeModel tree, MatrixParameter mu, Parameter indicators, Parameter locationDrift) {
    public DriftedMuStatistic( TreeModel tree, MatrixParameter mu, Parameter indicators,  Parameter mu1Scale, Parameter mu2Scale, Parameter muMean) {
        
        this.treeModel = tree;
        this.mu = mu;
        this.indicators = indicators;
       // this.locationDrift = locationDrift;
        this.mu1ScaleParameter = mu1Scale;
        mu1ScaleParameter.addParameterListener(this);
        this.mu2ScaleParameter = mu2Scale;
        mu2ScaleParameter.addParameterListener(this);
        
        this.muMeanParameter = muMean;
        muMeanParameter.addParameterListener(this);
        
        mu.addParameterListener(this);
        indicators.addParameterListener(this);
        //locationDrift.addParameterListener(this);
    }
    


    public int getDimension() {
        return treeModel.getNodeCount()*2;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {

    	int curNode = dim/2;
    	double value = mu.getParameter(curNode).getParameterValue(dim % 2);
    	
    	//if((int) indicators.getParameterValue(curNode)  == 0){
    	//	value = 0;
    	//}
    	

	    	if(  dim % 2 == 0 ){
	    		//value = value * locationDrift.getParameterValue(0) ;
	    		value =   value * mu1ScaleParameter.getParameterValue(0);
	    	}
	    	else{
	    		value = value * mu2ScaleParameter.getParameterValue(0);
			
	    	}
    
    	
          return (  value );

    }

    
    
     
    
    
    public String getDimensionName(int dim) {
    	String name = "mu_" +  ((dim/2) )  + "-" + ((dim %2 ) +1 );
        return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String MU_STRING = "mu";
        public final static String INDICATORS_STRING = "indicators";
       // public final static String LOCATION_DRIFT = "locationDrift";
        public final static String MU1_SCALE_PARAMETER = "mu1Scale";
        public final static String MU2_SCALE_PARAMETER = "mu2Scale";
        public final static String MU_MEAN_PARAMETER = "muMean";
        
        public String getParserName() {
            return DRIFTED_MU_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            MatrixParameter muParam = (MatrixParameter) xo.getElementFirstChild(MU_STRING);
            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS_STRING);
            //Parameter locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            Parameter mu1Scale = null;
            if (xo.hasChildNamed(MU1_SCALE_PARAMETER)) {
            	mu1Scale = (Parameter) xo.getElementFirstChild(MU1_SCALE_PARAMETER);
            }
            
            Parameter mu2Scale = null;
            if (xo.hasChildNamed(MU2_SCALE_PARAMETER)) {
            	mu2Scale = (Parameter) xo.getElementFirstChild(MU2_SCALE_PARAMETER);
            }  
            
            Parameter muMean = null;
            if(xo.hasChildNamed(MU_MEAN_PARAMETER)){
            	muMean = (Parameter) xo.getElementFirstChild(MU_MEAN_PARAMETER);
            }
            
            //return new DriftedMuStatistic( treeModel, muParam, indicators, locationDrift);
            return new DriftedMuStatistic( treeModel, muParam, indicators, mu1Scale, mu2Scale, muMean);
            

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return DriftedMuStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(MU_STRING, Parameter.class),
            new ElementRule(INDICATORS_STRING, Parameter.class),
            //new ElementRule(LOCATION_DRIFT, Parameter.class)
            new ElementRule(MU1_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the first dimension of mu"),
            new ElementRule(MU2_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the second dimension of mu"), 
            new ElementRule(MU_MEAN_PARAMETER, Parameter.class)
        };
    };

    

}
