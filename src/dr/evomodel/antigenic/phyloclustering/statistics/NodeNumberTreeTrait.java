package dr.evomodel.antigenic.phyloclustering.statistics;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class NodeNumberTreeTrait implements TreeTraitProvider {
	
    protected Helper treeTraits = new Helper();
    private TreeModel treeModel;
    public static final String NODE_NUMBER_TREETRAIT = "nodeNumberTreeTrait";

	
	 public NodeNumberTreeTrait(TreeModel treeModel_in){
		 
		 this.treeModel = treeModel_in;
		 
	        treeTraits.addTrait(new TreeTrait.IA() {
	        	
	            public String getTraitName() {
	            	return "node";
	            }

	            public String getTraitString(Tree tree, NodeRef node) {
	            	
	            	if(tree != treeModel){
	            		System.out.println("Something is wrong. Why is tree not equal to treeModel?");
	            		System.exit(0);
	            	}

	            	String nodeString = node.getNumber() + "";  //to get the node numbering of the tree.
	            	return nodeString;
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
	            	System.out.println("getTrait ran. Not expected. Quit now");
	            	System.exit(0);
	            	return null;
	            }


	        });

		 
	 }
	 
	 

	    public TreeTrait[] getTreeTraits() {
	        return treeTraits.getTreeTraits();
	    }
	
	
	    public TreeTrait getTreeTrait(String key) {
	    	System.out.println("not expected to run getTreeTrait. Quit now");
	    	System.exit(0);
	        return treeTraits.getTreeTrait(key);
	    }
	 

	    

	    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


	        public String getParserName() {
	            return NODE_NUMBER_TREETRAIT;
	        }

	        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

	        	 return new NodeNumberTreeTrait( treeModel);

	        }

	        //************************************************************************
	        // AbstractXMLObjectParser implementation
	        //************************************************************************

	        public String getParserDescription() {
	            return "Display node number in the tree.";
	        }

	        public Class getReturnType() {
	            return NodeNumberTreeTrait.class;
	        }

	        public XMLSyntaxRule[] getSyntaxRules() {
	            return rules;
	        }

	        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            new ElementRule(TreeModel.class),
	        };
	    };
	
	
}
