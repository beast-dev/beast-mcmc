
package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.evomodel.antigenic.phyloclustering.TreeClusteringVirusesPrior;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */


public class DriverCountStatistic extends Statistic.Abstract implements VariableListener {

	
    public static final String DRIVERCOUNT_STATISTIC = "driverCountStatistic";

    private TreeModel treeModel;
    private TreeClusteringVirusesPrior clusterPrior;


    
    public DriverCountStatistic(TreeModel tree, TreeClusteringVirusesPrior clusterPrior_in) {
        this.treeModel = tree;
		 this.clusterPrior = clusterPrior_in;
    }
    


    public int getDimension() {
        int[] causalCount = clusterPrior.getCausalCount();
		int numdata = causalCount.length *2;
        return numdata;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {
       	
       int[] causalCount = clusterPrior.getCausalCount();
       int[] nonCausalCount = clusterPrior.getNonCausalCount();
             
       int index =  dim/2;
       double value = -1;
   		if(dim%2==0){
   			value = causalCount[index];
   		}
   		else{
			value = nonCausalCount[index];
		}
        //System.out.println("dim=" + dim +  " dim%2=" + dim%2 + " and index = " + index + " value = " + value);

       return ( value);
    }


    
    public String getDimensionName(int dim) {
    	String name = "";
    	if(dim%2==0){
    		name += "C";
    	}
    	else{
    		name += "N";
    	}
    	name += ""+  (dim/2 +1);
    	return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {



        public String getParserName() {
            return DRIVERCOUNT_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            TreeClusteringVirusesPrior clusterPrior = (TreeClusteringVirusesPrior) xo.getChild(TreeClusteringVirusesPrior.class);
 
            return new DriverCountStatistic( treeModel, clusterPrior);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return ".";
        }

        public Class getReturnType() {
            return DriverCountStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            new ElementRule(TreeModel.class),
	            new ElementRule(TreeClusteringVirusesPrior.class),
        };
    };

    

}
