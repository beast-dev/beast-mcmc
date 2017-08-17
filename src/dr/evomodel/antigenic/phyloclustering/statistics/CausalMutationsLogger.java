package dr.evomodel.antigenic.phyloclustering.statistics;

import java.util.Iterator;
import java.util.LinkedList;

import dr.evomodel.antigenic.phyloclustering.TreeClusteringVirusesPrior;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class CausalMutationsLogger implements Loggable {

    /**
     * @return the log columns.
     */
    public LogColumn[] getColumns() {
    	 	
    	
    	int numDimension = treeModel.getNodeCount();
        LogColumn[] columns = new LogColumn[numDimension];
        for (int i = 0; i < numDimension; i++) {
            //columns[i] = new StatisticColumn(getDimensionName(i), i);
        	String traitName = "" + (i) + ":";
	    	LinkedList<Integer>[] mutationList= clusterPrior.getMutationList();
	    	if(mutationList[i] != null){
		    	Iterator itr = mutationList[i].iterator();
		    	int count = 0;
		    	while(itr.hasNext()){
		    		if(count > 0){
		    			traitName +=",";
		    		}
		    		int curMutation = ((Integer) itr.next()).intValue();
		    		traitName += curMutation;
		    		count++;
		    	}
	    	}
        	
        	
        	final int curNode = i;
        	columns[i] = new LogColumn.Abstract(traitName){
        		@Override
        		protected String getFormattedValue(){
        			//return "AAA";
        	    	LinkedList<Integer>[] causalList= clusterPrior.getCausalList();
        			return parseCausalList(causalList[curNode]);
        		}
        	};
        }
        return columns;
    }
    
    
    
    public String parseCausalList(LinkedList<Integer> causalMutations){
    	String stateList = "";
    	if(causalMutations != null){
	    	Iterator itr = causalMutations.iterator();
	    	while(itr.hasNext()){
	    		int curState = ((Integer) itr.next()).intValue();
	    		stateList += curState;
	    	}
	    	stateList = "s" + new StringBuilder(stateList).reverse().toString(); //need to reverse printing the states to be correct.
    	}
    	else{
    		stateList = "s";
    	}
    	 return(stateList);
    }
    
    public static final String PARSER_NAME = "causalMutationsLogger";
    private TreeModel treeModel;
    private TreeClusteringVirusesPrior clusterPrior;
    
    public CausalMutationsLogger(TreeModel tree, TreeClusteringVirusesPrior clusterPrior_in) {
        this.treeModel = tree;
		 this.clusterPrior = clusterPrior_in;
    }
    
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {



        public String getParserName() {
            return PARSER_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            TreeClusteringVirusesPrior clusterPrior = (TreeClusteringVirusesPrior) xo.getChild(TreeClusteringVirusesPrior.class);
 
            return new CausalMutationsLogger( treeModel, clusterPrior);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return ".";
        }

        public Class getReturnType() {
            return CausalMutationsLogger.class;
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
