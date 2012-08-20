package dr.evomodel.speciation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


import test.dr.evomodel.speciation.AlloppSpeciesNetworkModelTEST;
import test.dr.evomodel.speciation.AlloppSpeciesNetworkModelTEST.LogLhoodGTreeInNetworkTEST;

import jebl.util.FixedBitSet;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.MutableTreeListener;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeLogger;
import dr.math.MathUtils;
import dr.util.AlloppMisc;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;


/**
 * 
 * A tree for a single ploidy level in an allopolyploid network.
 * 
 * @author Graham Jones
 *         Date: 01/05/2011
 */



/*
 * class AlloppLeggedTree
 * 
 * This is a `tree with legs', which is a homoploid 
 * species tree which is attached to a tree of lower ploidy
 * via its legs, as part of a AlloppSpeciesNetworkModel. 
 * 
 * 2011-05-19 I use this for diploid tree, although only the tree,
 * no legs, is used then. 
 * 
 * tree is a SimpleTree: its nodes contain times, and taxa at tips.
 * 
 * For tetraploids, legs[] has length one or two. If two legs, they specify
 * the branches (which may be the same branch twice at different
 * times) where the tree joins a lower ploidy tree.
 * 
 * If only one leg, it means that the polyploid
 * tree arose from two extinct species. In this
 * case the field splitheight is used for the MRCA of these,
 * and the single leg joins the MRCA to the lower ploidy tree. 
 * A special case of the one leg case is where there is
 * no diploid tree. Here the leg has no foot.
 * 
 * hybridheight is the time of hybridization, `where the legs 
 * join together', which is earlier then the root node of the 
 * homoploid tree.
 * 
 */




public class AlloppLeggedTree implements MutableTree, TreeLogger.LogUpon  {
	
    private  SimpleTree tree;
    private  Leg[] legs;
    private double splitheight;
    private double hybridheight;

    
    
    private class Leg {
    	/*
    	 *  footUnion is the node in a lower ploidy tree Y whose branch 
    	 *  contains the foot leading to this tree. footUnion specifies
    	 *  the clade (of species) in Y at the node.
    	 *  
    	 *  grjtodo tetraonly. With hexaploids, etc, may need to identify Y?
    	 *  Or clade is enough?
    	 *  
    	 *  height is the time during the branch of the foot. 
    	 *  
    	 *  grjtodo? footID identifies the foot - heights might be identical 
    	 */
    	public FixedBitSet footUnion;
    	public double height;
    	
    	/**
    	 * clone constructor
    	 */
    	public Leg(Leg leg) {
    		this.height = leg.height;
    		this.footUnion = new FixedBitSet(leg.footUnion);
    	}
    	
    	// Partial constructor. The leg dangles, unattached
    	public Leg(double height) {
    		this.footUnion = new FixedBitSet(0);
    		this.height = height;
    	}
    }

    
    
    /*
     * Constructor makes a random starting (homoploid) tree with legs.
     */
    public AlloppLeggedTree(Taxon[] taxa, AlloppSpeciesNetworkModel.LegType legtype) {
    	SimpleNode root;
    	int nTaxa = taxa.length;
    	int nNodes = 2 * nTaxa - 1;
    	SimpleNode[] nodes = new SimpleNode[nNodes];
    	for (int n = 0; n < nNodes; n++) {
    		nodes[n] = new SimpleNode();
    	}    
    	ArrayList<Integer> tojoin = new ArrayList<Integer>(nTaxa);
    	for (int n = 0; n < nTaxa; n++) {
    		nodes[n].setTaxon(taxa[n]);
    		nodes[n].setHeight(0.0);
    		tojoin.add(n);
    	}
    	double rate = 1.0;
    	double treeheight = 0.0;
    	for (int i = 0; i < nTaxa-1; i++) {
    		int numtojoin = tojoin.size();
    		int j = MathUtils.nextInt(numtojoin);
    		Integer child0 = tojoin.get(j);
    		tojoin.remove(j);
    		int k = MathUtils.nextInt(numtojoin-1);
    		Integer child1 = tojoin.get(k);
    		tojoin.remove(k);
    		nodes[nTaxa+i].addChild(nodes[child0]);
    		nodes[nTaxa+i].addChild(nodes[child1]);
    		nodes[nTaxa+i].setHeight(treeheight + randomnodeheight(numtojoin*rate));
    		treeheight = nodes[nTaxa+i].getHeight();
    		tojoin.add(nTaxa+i);
    	}
    	hybridheight = -1.0;
    	splitheight = -1.0;
    	switch (legtype) {
    	case NODIPLOIDS:
    		hybridheight = treeheight + randomnodeheight(rate);
    		splitheight =  hybridheight + randomsplitheight(rate);
    		legs = new Leg[1];
    		legs[0] = new Leg(-1.0); // no foot
    		break;
    	case JOINED:
    		hybridheight = treeheight + randomnodeheight(rate);
    		splitheight =  hybridheight + randomsplitheight(rate);
    		legs = new Leg[1];
    		legs[0] = new Leg(splitheight + randomnodeheight(rate));
    		break;
    	case ONEBRANCH: case TWOBRANCH:
    		hybridheight = treeheight + randomnodeheight(rate);
    		legs = new Leg[2];
    		legs[0] = new Leg(hybridheight + randomnodeheight(rate));
    		legs[1] = new Leg(legs[0].height + randomnodeheight(rate));
    		break;
    	case NONE:
    		legs = new Leg[0];
    		break;
    	default:
    		assert false;
    		legs = new Leg[0];
    	}
    	root = nodes[nodes.length - 1];
    	tree = new SimpleTree(root);
    	tree.setUnits(Units.Type.SUBSTITUTIONS); 
    }
    
   
    /** 
     * clone constructor 
     */
    public AlloppLeggedTree(AlloppLeggedTree tree) {
    	this.tree = new SimpleTree(tree);
    	this.hybridheight = tree.hybridheight;
    	this.splitheight = tree.splitheight;
    	this.legs = new Leg[tree.legs.length];
    	for (int i = 0; i < legs.length; i++) {
    		legs[i] = new Leg(tree.legs[i]);
    	}
    }
     
    
     /*
      * Constructor for testing. 
      */
     public AlloppLeggedTree(Taxon[] taxa, AlloppSpeciesNetworkModelTEST.NetworkToMulLabTreeTEST nmltTEST,
    		                       AlloppSpeciesNetworkModel.LegType legtype, double addheight) {
         int nTaxa = taxa.length;
         int nNodes = 2 * nTaxa - 1;
         SimpleNode[] nodes = new SimpleNode[nNodes];
         for (int n = 0; n < nNodes; n++) {
             nodes[n] = new SimpleNode();
         }
         SimpleNode root;
         nodes[0].setTaxon(taxa[0]);
         if (nTaxa >= 2) {
         nodes[1].setTaxon(taxa[1]);
         nodes[2].setHeight(addheight + 1.0);
         nodes[2].addChild(nodes[0]);
         nodes[2].addChild(nodes[1]);
         }
         if (nTaxa == 3) {
             nodes[3].setTaxon(taxa[2]);
             nodes[4].setHeight(addheight + nodes[2].getHeight() + 1.0);
             nodes[4].addChild(nodes[2]);
             nodes[4].addChild(nodes[3]);            
         }
         root = nodes[nodes.length - 1];
         tree = new SimpleTree(root);
         tree.setUnits(Units.Type.YEARS);
         
         double rootheight = root.getHeight();
         switch (legtype) {
         case NONE:
        	 legs = new Leg[0];
        	 splitheight = -1.0;
        	 break;
         case TWOBRANCH: case ONEBRANCH:
        	 legs = new Leg[2];
        	 legs[0] = new Leg(rootheight+1.0);
        	 legs[1] = new Leg(rootheight+2.0);
        	 splitheight = -1.0;
        	 break;
         case JOINED:
        	 legs = new Leg[1];
        	 legs[0] = new Leg(rootheight+2.0);
        	 splitheight = rootheight+1.0;
        	 break;
         default:
        	 assert false;
        	 legs = new Leg[0];
         }
     }       
     
     
     /*
      * Constructor for testing
      */
     public AlloppLeggedTree(Taxon[] taxa, AlloppSpeciesNetworkModel.LegType legtype,
			LogLhoodGTreeInNetworkTEST llgtnTEST) {
    	 assert legtype == AlloppSpeciesNetworkModel.LegType.NODIPLOIDS;
    	 int nTaxa = taxa.length;
    	 assert nTaxa <= 3;
    	 
    	 if (llgtnTEST.heights.length > 0) { // testing cases
    		 assert (nTaxa == 2 && llgtnTEST.heights.length == 3)  ||  (nTaxa == 3 && llgtnTEST.heights.length == 4);
    		 hybridheight = llgtnTEST.heights[nTaxa-1];
    		 splitheight = llgtnTEST.heights[nTaxa];
    	 }
    	 
         int nNodes = 2 * nTaxa - 1;
         SimpleNode[] nodes = new SimpleNode[nNodes];
         for (int n = 0; n < nNodes; n++) {
             nodes[n] = new SimpleNode();
         }
         SimpleNode root;
         nodes[0].setTaxon(taxa[0]);
         if (nTaxa >= 2) {
         nodes[1].setTaxon(taxa[1]);
         nodes[2].setHeight(llgtnTEST.heights[0]);
         nodes[2].addChild(nodes[0]);
         nodes[2].addChild(nodes[1]);
         }
         if (nTaxa == 3) {
             nodes[3].setTaxon(taxa[2]);
             nodes[4].setHeight(llgtnTEST.heights[1]);
             nodes[4].addChild(nodes[2]);
             nodes[4].addChild(nodes[3]);            
         }
         root = nodes[nodes.length - 1];
         tree = new SimpleTree(root);
         tree.setUnits(Units.Type.SUBSTITUTIONS);
         
    	 legs = new Leg[1];
    	 legs[0] = new Leg(-1.0); // no foot

	}

   
     
     
     public int scaleAllHeights(double scale) {
    	beginTreeEdit();
     	int nNodes = tree.getNodeCount();
     	for (int n = 0; n < nNodes; n++) {
     		NodeRef node = tree.getNode(n);
     		tree.setNodeHeight(node, scale*tree.getNodeHeight(node));
     	}
     	int count = nNodes;
     	if (hybridheight >= 0.0) { hybridheight *= scale;  count++; }
     	if (splitheight >= 0.0) { splitheight *= scale;  count++; }
     	for (int i = 0; i < legs.length; i++) {
     		if (legs[i].height >= 0.0) { legs[i].height *= scale;  count++; }
     	}
     	return count;
     }
     

     public double[] getInternalHeights() {
    	 int n = tree.getInternalNodeCount();
    	 double heights[] = new double[n];
    	 for (int i = 0; i < n; i++) {
    		 NodeRef node = tree.getInternalNode(i);
    		 heights[i] = tree.getNodeHeight(node);
    	 }    	  
    	 return heights;
     }   
     
     public double getRootHeight() {
    	 return tree.getRootHeight();
     }     

     public double getSplitHeight() {
    	 return splitheight;
     }


     public double getHybridHeight() {
    	 return hybridheight;
     }

     // 2011-08-31 used only for setting up cases for testing 
     public double getMaxFootHeight() {
    	 double max = tree.getRootHeight();
    	 for (int i = 0; i < legs.length; i++) {
    		 max = Math.max(max, legs[i].height);
    	 }
    	 return max;
     }
     
     
     // Returns the 'bottom' of a tree with legs. Works for 0,1,2 legs.
     // Used for making starting tree
     public double getMaxHeight() {
    	 double max = tree.getRootHeight();
    	 for (int i = 0; i < legs.length; i++) {
    		 max = Math.max(max, legs[i].height);
    	 }
    	 max = Math.max(max, hybridheight);
    	 max = Math.max(max, splitheight);
    	 return max;
     }


    
     
	public void setFootUnion(int leg, FixedBitSet footUnion) {
		legs[leg].footUnion = footUnion;
	}
	
	

	public FixedBitSet getFootUnion(int leg) {
		return legs[leg].footUnion;
	}
 
	public int getNumberOfLegs() {
		return legs.length;
	}

	public double getFootHeight(int leg) {
		return legs[leg].height;
	}

	
	public void setHybridHeight(double newh) {
		assert newh >= tree.getRootHeight();
		if (legs.length == 1) {
			assert newh <= splitheight;
		} else {
			assert newh <= legs[0].height;
			assert newh <= legs[1].height;
		}
		hybridheight = newh;
	}



	
	public void setSplitHeight(double news) {
		assert legs.length == 1;
		assert legs[0].height < 0.0  ||  news <= legs[0].height;
		assert news >= hybridheight;
		splitheight = news;
	}

	
	
	
    
    // Moves 

    public void moveSplitOrLeg(double dirooth) {
    	if (getNumberOfLegs() == 1) {
    		if (MathUtils.nextBoolean() == true) {
    			moveSplitHeight(dirooth);
    		} else {
    			moveLegHeight(0, dirooth);
    		}
    	} else {
    		if (MathUtils.nextBoolean() == true) {
    			moveLegHeight(0, dirooth);
    		} else {
    			moveLegHeight(0, dirooth);
    		}    		 
    	}

    }


    // tetraonly twodiploidsonly     
    public void moveMostRecentLegHeight() {
    	if (MathUtils.nextBoolean() == true) {
    		// small move
    		if (getNumberOfLegs() == 1) {
    			double oldh = getSplitHeight();
    			double minh = getHybridHeight();
    			double maxh = legs[0].height;
    			setSplitHeight(AlloppMisc.uniformInRange(oldh, minh, maxh, 0.1));
    		} else {
    			assert getNumberOfLegs() == 2;
    			if (legs[0].height < legs[1].height) {
    				legs[0].height = AlloppMisc.uniformInRange(legs[0].height, getHybridHeight(), legs[1].height, 0.1);
    			} else {
    				legs[1].height = AlloppMisc.uniformInRange(legs[1].height, getHybridHeight(), legs[0].height, 0.1); 
    			}
    		}      		 
    	} else {
    		// big move
    		if (getNumberOfLegs() == 1) {
    			setSplitHeight(MathUtils.uniform(getHybridHeight(), legs[0].height));
    		} else {
    			assert getNumberOfLegs() == 2;
    			if (legs[0].height < legs[1].height) {
    				legs[0].height = MathUtils.uniform(getHybridHeight(), legs[1].height);
    			} else {
    				legs[1].height = MathUtils.uniform(getHybridHeight(), legs[0].height);
    			}
    		}      		 
    	}

    }

    
    // tetraonly twodiploidsonly
    public void moveMostAncientLegHeight(double dirooth) {
    	if (MathUtils.nextBoolean() == true) {
    		// small move
    		if (getNumberOfLegs() == 1) {
    			legs[0].height = AlloppMisc.uniformInRange(legs[0].height, getSplitHeight(), dirooth, 0.1);
    		} else {
    			assert getNumberOfLegs() == 2;
    			if (legs[0].height < legs[1].height) {
    				legs[1].height = AlloppMisc.uniformInRange(legs[1].height, legs[0].height, dirooth, 0.1);
    			} else {
    				legs[0].height = AlloppMisc.uniformInRange(legs[0].height, legs[1].height, dirooth, 0.1);
    			}
    		}     		 
    	} else {
    		// big move
    		if (getNumberOfLegs() == 1) {
    			legs[0].height = MathUtils.uniform(getSplitHeight(), dirooth);
    		} else {
    			assert getNumberOfLegs() == 2;
    			if (legs[0].height < legs[1].height) {
    				legs[1].height = MathUtils.uniform(legs[0].height, dirooth);
    			} else {
    				legs[0].height = MathUtils.uniform(legs[1].height, dirooth);
    			}
    		}    		 
    	}

    }


    // grjtodo tetraonly twodiploidsonly
     public void moveLegTopology(FixedBitSet dip0, FixedBitSet dip1) {
    	 int rnd = MathUtils.nextInt(6);
    	 double t0;
    	 double t1;
    	 if (getNumberOfLegs() == 1) {
    		 t0 = getSplitHeight();
    		 t1 = legs[0].height;
    	 } else {
    		 t0 = Math.min(legs[0].height, legs[1].height);
    		 t1 = Math.max(legs[0].height, legs[1].height);
    	 }
    	 switch (rnd) {
    	 case 0:  case 1:  case 2:  case 3:
    		 // TWOBRANCH, ONEBRANCH
    		 legs = new Leg[2];
    		 if (MathUtils.nextBoolean()) {
    		   legs[0] = new Leg(t0);
    		   legs[1] = new Leg(t1);
    		 } else {
    		   legs[0] = new Leg(t1);
    		   legs[1] = new Leg(t0);
    		 }
    			 
    		 break;
    	 case 4: case 5:
    		 // JOINED
    		 legs = new Leg[1];
    		 legs[0] = new Leg(t1);
    		 setSplitHeight(t0);
    		 break;
    	 }
    	 switch (rnd) {
    	 case 0:
    		 legs[0].footUnion = dip0;
    		 legs[1].footUnion = dip1;
    		 break;
    	 case 1:
    		 legs[0].footUnion = dip1;
    		 legs[1].footUnion = dip0;
    		 break;
    	 case 2:
    		 legs[0].footUnion = dip0;
    		 legs[1].footUnion = dip0;
    		 break;
    	 case 3:
    		 legs[0].footUnion = dip1;
    		 legs[1].footUnion = dip1;
    		 break;
    	 case 4:
    		 legs[0].footUnion = dip0;
    		 break;
    	 case 5:
    		 legs[0].footUnion = dip1;
    		 break;
    	 }    	 
     }
	
	
	
	
	
	
	
	
	
	
	

	// TreeLogger.LogUpon
	
	public boolean logNow(long state) {
       //		grjtodo 

		return false;
	}


	
	
	
	
/*
 * 
 * **************************************************************
 * 
 * 	
 */
	
	
	
	private void moveLegHeight(int lg, double dirooth) {
		double minh = getHybridHeight();
		if (getNumberOfLegs() == 1) {
			minh = getSplitHeight();
		}
		if (MathUtils.nextBoolean() == true) {
			// small move
			legs[lg].height = AlloppMisc.uniformInRange(legs[lg].height, minh, dirooth, 0.1);
		} else {
			// big move
			legs[lg].height = MathUtils.uniform(minh, dirooth);
		}
	}


	private void moveSplitHeight(double dirooth) {
		double oldh = getSplitHeight();
		double minh = getHybridHeight();
		double maxh = legs[0].height;
		if (MathUtils.nextBoolean() == true) {
			// small move
			setSplitHeight(AlloppMisc.uniformInRange(oldh, minh, maxh, 0.1));
		} else {
			// big move
			setSplitHeight(MathUtils.uniform(getHybridHeight(), legs[0].height));
		}
	}
	
	
	
    private double randomnodeheight(double rate) {
    	return MathUtils.nextExponential(rate) + 1e-6/rate;
    	// 1e-6/rate to avoid very tiny heights
    }
    
    
    private double randomsplitheight(double rate) {
    	return randomnodeheight(rate) + randomnodeheight(rate) + randomnodeheight(rate);
     }	
	
    

//                      Tree   
    
	
	public NodeRef getRoot() {
		return tree.getRoot();
	}

	
	public int getNodeCount() {
		return tree.getNodeCount();    
	}

	
	public NodeRef getNode(int i) {
		return tree.getNode(i);    
	}

	
	public NodeRef getInternalNode(int i) {
		return tree.getInternalNode(i);    
	}

	
	public NodeRef getExternalNode(int i) {
		return tree.getExternalNode(i);    
	}

	
	public int getExternalNodeCount() {
		return tree.getExternalNodeCount();    
	}

	
	public int getInternalNodeCount() {
		return tree.getInternalNodeCount();    
	}

	
	public Taxon getNodeTaxon(NodeRef node) {
		return tree.getNodeTaxon(node);    
	}

	
	public boolean hasNodeHeights() {
		return tree.hasNodeHeights();    
	}

	
	public double getNodeHeight(NodeRef node) {
		return tree.getNodeHeight(node);    
	}

	
	public boolean hasBranchLengths() {
		return tree.hasBranchLengths();    
	}

	
	public double getBranchLength(NodeRef node) {
		return tree.getBranchLength(node);    
	}

	
	public double getNodeRate(NodeRef node) {
		return tree.getNodeRate(node);    
	}

	
	public Object getNodeAttribute(NodeRef node, String name) {
		return tree.getNodeAttribute(node, name);    
	}

	
	public Iterator getNodeAttributeNames(NodeRef node) {
		return tree.getNodeAttributeNames(node);    
	}

	
	public boolean isExternal(NodeRef node) {
		return tree.isExternal(node);    
	}

	
	public boolean isRoot(NodeRef node) {
		return tree.isRoot(node);    
	}

	
	public int getChildCount(NodeRef node) {
		return tree.getChildCount(node);    
	}

	
	public NodeRef getChild(NodeRef node, int j) {
		return tree.getChild(node, j);    
	}

	
	public NodeRef getParent(NodeRef node) {
		return tree.getParent(node);    
	}

	
	public Tree getCopy() {
		return tree.getCopy();    
	}

	
	public int getTaxonCount() {
		return tree.getTaxonCount();    
	}

	
	public Taxon getTaxon(int taxonIndex) {
		return tree.getTaxon(taxonIndex);    
	}

	
	public String getTaxonId(int taxonIndex) {
		return tree.getTaxonId(taxonIndex);    
	}

	
	public int getTaxonIndex(String id) {
		return tree.getTaxonIndex(id);    
	}

	
	public int getTaxonIndex(Taxon taxon) {
		return tree.getTaxonIndex(taxon);    
	}

	
	public List<Taxon> asList() {
		return tree.asList();    
	}

	
	public Object getTaxonAttribute(int taxonIndex, String name) {
		return tree.getTaxonAttribute(taxonIndex, name);    
	}

	
	public Iterator<Taxon> iterator() {
		return tree.iterator();    
	}

	
	public Type getUnits() {
		return tree.getUnits();    
	}

	
	public void setUnits(Type units) {
		tree.setUnits(units);    
	}

	
	public void setAttribute(String name, Object value) {
		tree.setAttribute(name, value);    
	}

	
	public Object getAttribute(String name) {
		return tree.getAttribute(name);    
	}

	
	public Iterator<String> getAttributeNames() {
		return tree.getAttributeNames();    
	}

	
	public int addTaxon(Taxon taxon) {
		return tree.addTaxon(taxon);    
	}

	
	public boolean removeTaxon(Taxon taxon) {
		return tree.removeTaxon(taxon);    
	}

	
	public void setTaxonId(int taxonIndex, String id) {
		tree.setTaxonId(taxonIndex, id);
	}

	
	public void setTaxonAttribute(int taxonIndex, String name, Object value) {
		tree.setTaxonAttribute(taxonIndex, name, value);
	}

	
	public void addMutableTaxonListListener(MutableTaxonListListener listener) {
		tree.addMutableTaxonListListener(listener);
	}

	
	
	//          MutableTree which extends Tree, MutableTaxonList


		
	public boolean beginTreeEdit() {
		return tree.beginTreeEdit();
	}

	public void endTreeEdit() {
		tree.endTreeEdit();
	}

	
	public void addChild(NodeRef parent, NodeRef child) {
		tree.addChild(parent, child);
	}

	public void removeChild(NodeRef parent, NodeRef child) {
		tree.removeChild(parent, child);
	}

	public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
		tree.replaceChild(node, child, newChild);
	}

	public void setRoot(NodeRef root) {
		tree.setRoot(root);
	}

	
	public void setNodeHeight(NodeRef node, double height) {
		tree.setNodeHeight(node, height);
	}

	
	public void setNodeRate(NodeRef node, double rate) {
		tree.setNodeRate(node, rate);
	}

	
	public void setBranchLength(NodeRef node, double length) {
		tree.setBranchLength(node, length);
	}

	
	public void setNodeAttribute(NodeRef node, String name, Object value) {
		tree.setNodeAttribute(node, name, value);
	}

	
	public void addMutableTreeListener(MutableTreeListener listener) {
		tree.addMutableTreeListener(listener);
	}

	
	
	
	
	// Identifiable

	public String getId() {
		return tree.getId();
	}



	
	public void setId(String id) {
		tree.setId(id);
		
	}


}
