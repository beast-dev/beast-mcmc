package dr.evomodel.antigenic.phyloclustering.statistics;

import dr.evolution.datatype.Nucleotides;

import java.util.Iterator;
import java.util.LinkedList;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.antigenic.phyloclustering.TreeClusteringVirusesPrior;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class MutationsTreeTrait implements TreeTraitProvider {
	
    protected Helper treeTraits = new Helper();
    private TreeModel treeModel;
    private TreeClusteringVirusesPrior clusterPrior;
    public static final String MUTATIONS_TREETRAIT = "mutationsTreeTrait";

    
    private String[] mutationString;
    
    private LinkedList<Integer>[] mutationList;
    private LinkedList<Integer>[] causalList;
    
	
	 public MutationsTreeTrait(TreeModel treeModel_in, TreeClusteringVirusesPrior clusterPrior_in){
		 
		 this.treeModel = treeModel_in;
		 this.clusterPrior = clusterPrior_in;
		 

		 		 
	        //alignment.setDataType(siteModel.getSubstitutionModel().getDataType());

		 	int numNodes = treeModel.getNodeCount();
	        // Get sequences
	        String[] sequence = new String[numNodes];
	        
	     // Universal
    		String GENETIC_CODE_TABLES ="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

    		int numCodons = clusterPrior.getNumSites();
    		 //int numCodons = 330;
	        for(int curIndex = 0; curIndex < numNodes; curIndex ++){
	    		String ns =  (String) treeModel.getNodeAttribute( treeModel.getNode(curIndex), "states");
	
	    		ns = ns.substring(clusterPrior.getStartBase(), clusterPrior.getEndBase() );
	    		//ns = ns.substring(3+27, ns.length() - 1);
	    		//System.out.println(ns);
	    		
	    		
	    		//numCodons = ns.length()/3;  // or do I care about only 330?

	    		//System.out.println(numCodons);
	    		String codonSequence = "";
	    		for(int codon=0; codon< numCodons; codon++){
	    			
	    			int nuc1 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3)];
	    			int nuc2 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3+1)];
	    			int nuc3 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3+2)];
	    			
	    			int canonicalState = (nuc1 * 16) + (nuc2 * 4) + nuc3;
	    			
	    			codonSequence = codonSequence + GENETIC_CODE_TABLES.charAt(canonicalState);
	    		}
				//System.out.println(codonSequence);
	            sequence[curIndex] = codonSequence;
	    		
	        }


	        mutationString = new String[treeModel.getNodeCount()];

			NodeRef cNode = treeModel.getRoot();
		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
		    
		    visitlist.add(cNode);
		    
		    int countProcessed=0;
		    while(visitlist.size() > 0){
		    	countProcessed++;
		    	//assign value to the current node...
		    	if(treeModel.getParent(cNode) == null){  //this means it is a root node
		    		//visiting the root
		    		//System.out.println(cNode.getNumber() + ":\t" + "root");
		    	}
		    	else{
		    		//visiting
		    		//System.out.print(cNode.getNumber() + ":\t");

		    		//String listMutations = "\"";
		    		mutationString[cNode.getNumber()]  = "\"";
		    		String nodeState =  sequence[cNode.getNumber()];
		    		String parentState =  sequence[treeModel.getParent(cNode).getNumber()];
		    		           
		    		int count = 0;
		    		for(int i=0; i < numCodons; i++){
		    			if(nodeState.charAt(i) != parentState.charAt(i)){
		    				count++;
		    				if(count>1){
		    					//System.out.print(",");
		    					mutationString[cNode.getNumber()] =  mutationString[cNode.getNumber()] + ",";
		    				}
		    				//System.out.print(i+1);
		    				mutationString[cNode.getNumber()] =  mutationString[cNode.getNumber()] + (i+1);  //i+1 so mutation starts from 1 - 330
		    			}
		    			
		    			//store in linked list
		    		}
		    		//System.out.println("");
		    		mutationString[cNode.getNumber()]  = mutationString[cNode.getNumber()]  + "\"";
		    	}
		    	
				//System.out.println(cNode.getNumber() + "\t" +  treeModel.getNodeAttribute(cNode, "states") );

		    	
		    	//add all the children to the queue
	  			for(int childNum=0; childNum < treeModel.getChildCount(cNode); childNum++){
	  				NodeRef node= treeModel.getChild(cNode,childNum);
	  				visitlist.add(node);
	  	        }
	  			
		  			
		  		visitlist.pop(); //now that we have finished visiting this node, pops it out of the queue

	  			if(visitlist.size() > 0){
	  				cNode = visitlist.getFirst(); //set the new first node in the queue to visit
	  			}
	  			
				
		}

treeTraits.addTrait(new TreeTrait.IA() {
	
	        	
	            public String getTraitName() {          	
	            	return "mutations";
	            }

	            public String getTraitString(Tree tree, NodeRef node) {
	            	if(tree != treeModel){
	            		System.out.println("Something is wrong. Why is tree not equal to treeModel?");
	            		System.exit(0);
	            	}

	            	//String nodeString = mutationString[node.getNumber()];  //to get the node numbering of the tree.

	            	
	      		  	mutationList = clusterPrior.getMutationsPerNode();
	      		  	causalList = clusterPrior.getCausativeStatesPerNode();
	            	
	            	String nodeString = "\"";
	            			
        	    	if(mutationList[node.getNumber()] != null){
        		    	Iterator itr = mutationList[node.getNumber()].iterator();
        		    	Iterator itr2 = causalList[node.getNumber()].iterator();
        		    	int count = 0;
        		    	while(itr.hasNext()){
        		    		count++;
        		    		if(count>1){
        		    			nodeString += ",";
        		    		}
        		    		int curMutation = ((Integer) itr.next()).intValue();
        		    		nodeString += curMutation;
        		    		int curCausal = ((Integer) itr2.next()).intValue();
        		    		if(curCausal ==1){
        		    			nodeString += "*";
        		    		}
        		    	}
        	    	}
        	    	nodeString += "\"";
	            	
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
	            return MUTATIONS_TREETRAIT;
	        }

	        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
	            TreeClusteringVirusesPrior clusterPrior = (TreeClusteringVirusesPrior) xo.getChild(TreeClusteringVirusesPrior.class);
	        	 return new MutationsTreeTrait( treeModel, clusterPrior);

	        }

	        //************************************************************************
	        // AbstractXMLObjectParser implementation
	        //************************************************************************

	        public String getParserDescription() {
	            return "Display node number in the tree.";
	        }

	        public Class getReturnType() {
	            return MutationsTreeTrait.class;
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
