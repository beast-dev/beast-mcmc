package dr.evomodel.antigenic.phyloclustering.statistics;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ClusterLabelsVirusesTreeTrait implements TreeTraitProvider {
	
    protected Helper treeTraits = new Helper();
    private TreeModel treeModel;
    public static final String CLUSTERLABELSTREETRAIT = "clusterLabelsVirusesTreeTrait";
    public final static String CLUSTERLABELSTREENODE = "clusterLabelsTreeNodes";
    
    private Parameter clusterLabelsTreeNode;

	
	 public ClusterLabelsVirusesTreeTrait(TreeModel treeModel_in, Parameter clusterLabelsTreeNode_in){
		 
		 this.treeModel = treeModel_in;
		 this.clusterLabelsTreeNode = clusterLabelsTreeNode_in;
		 
	        treeTraits.addTrait(new TreeTrait.IA() {
	        	
	            public String getTraitName() {
	            	//System.out.println("print label");
	               // return tag;
	            	return "cluster";
	            }

	            public String getTraitString(Tree tree, NodeRef node) {
	            	
	            	if(tree != treeModel){
	            		System.out.println("Something is wrong. Why is tree not equal to treeModel?");
	            		System.exit(0);
	            	}

	            	//the problem is, I don't know how to only do the processing to get the clusterLabels
	            	//right before this routine is run... 
	            	//so I have to keep a parameter or variable somewhere to store the clusterLabels information of the tree nodes
	            	//whenever the cluster assignment is changed.. 
	            	//the states have to be precomputed so this only prints..
	            	
	            	//System.out.println("node=" + node.getNumber());
	            	
	            	String clusterLabelString = ((int) clusterLabelsTreeNode.getParameterValue(node.getNumber()) ) + "";
	            	//String clusterLabelString = node.getNumber() + "";  //to get the node numbering of the tree.
	            	
	            	//if(node.getNumber() ==0){
	            	//	//System.out.println("print ");
	            	//	for(int i=0; i < treeModel.getNodeCount(); i++){
	            	//		System.out.print(  clusterLabelsTreeNode.getParameterValue(i) + "\t");
	            	//	}
	            	//	System.out.println("");
	            	//}
	                //return formattedState(getStatesForNode(tree, node), dataType);
	            	return clusterLabelString;
	            }
	            
	            
	            public Intent getIntent() {
	            	//System.out.println("getIntent");
	                return Intent.NODE;
	            }

	            public Class getTraitClass() {
	            	System.out.println("getTraitClass ran. Not expected. Quit now");
	            	System.exit(0);
	                return int[].class;
	            }

	            
	            public int[] getTrait(Tree tree, NodeRef node) {
	              //  return getStatesForNode(tree, node);
	            	System.out.println("getTrait ran. Not expected. Quit now");
	            	System.exit(0);
	            	//int x[] = new int[10];
	            	return null;
	            }


	        });

		 
	 }
	 
	 

	    public TreeTrait[] getTreeTraits() {
	    	//System.out.println("hihi");
	        return treeTraits.getTreeTraits();
	    }
	
	
	    public TreeTrait getTreeTrait(String key) {
	    	System.out.println("not expected to run getTreeTrait. Quit now");
	    	System.exit(0);
	        return treeTraits.getTreeTrait(key);
	    }
	 

	    

	    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


	        public String getParserName() {
	            return CLUSTERLABELSTREETRAIT;
	        }

	        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

	            XMLObject cxo = xo.getChild(CLUSTERLABELSTREENODE);
                Parameter clusterLabelsTreeNode = (Parameter) cxo.getChild(Parameter.class);
           
	        	 return new ClusterLabelsVirusesTreeTrait( treeModel, clusterLabelsTreeNode);

	        }

	        //************************************************************************
	        // AbstractXMLObjectParser implementation
	        //************************************************************************

	        public String getParserDescription() {
	            return "Integrate ClusterLabels of viruses into the tree.";
	        }

	        public Class getReturnType() {
	            return ClusterLabelsVirusesTreeTrait.class;
	        }

	        public XMLSyntaxRule[] getSyntaxRules() {
	            return rules;
	        }

	        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            new ElementRule(TreeModel.class),
                new ElementRule(CLUSTERLABELSTREENODE, Parameter.class),

	        };
	    };
	
	
}
