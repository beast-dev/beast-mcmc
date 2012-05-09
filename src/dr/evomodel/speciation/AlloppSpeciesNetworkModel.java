package dr.evomodel.speciation;



import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeLogger;
import dr.evomodelxml.speciation.AlloppSpeciesNetworkModelParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.Scalable;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import jebl.util.FixedBitSet;
import java.util.*;
import java.util.logging.Logger;

import test.dr.evomodel.speciation.AlloppSpeciesNetworkModelTEST;


/**
 * 
 * Implements an allopolyploid species network as a collection of `trees with legs'.
 * 
 * @author Graham Jones
 *         Date: 19/04/2011
 */


/*
 * class AlloppSpeciesNetworkModel
 * 
 * Implements the species network as a collection of `trees with legs'.
 * and converts this representation into a multiply labelled 
 * binary tree.
 * 
 * General idea is that the network is easiest to change (eg detach
 * and re-attach tetraploid subtrees) while likelihood calculations
 * are easiest to do in the multiply labelled tree.
 * 
 * The individual `trees with legs' are implemented by AlloppLeggedTree's. 
 * The multiply labelled binary tree is implemented by AlloppMulLabTree.
 * 
 * *********************
 * 
 * apsp is a reference to the AlloppSpeciesBindings which knows how 
 * species are made of individuals and individuals are made of taxa,
 * and which contains the list of gene trees.
 * 
 * trees[][] represents the network as a set of homoploid trees
 * 
 * mullabtree represents the network as single tree with tips that
 * can be multiply labelled with species. 
 * 
 */
//  grjtodo JH's SpeciesTreeModel implements 
// MutableTree, TreeTraitProvider, TreeLogger.LogUpon, Scalable
// not clear how much of those sensible here.
// AlloppLeggedTree implements MutableTree, TreeLogger.LogUpon.
// Nothing so far does TreeTraitProvider.
public class AlloppSpeciesNetworkModel extends AbstractModel implements
		Scalable, Units, Citable, Tree, TreeLogger.LogUpon {

	private final AlloppSpeciesBindings apsp;
	private AlloppLeggedTree[][] trees;
	private AlloppLeggedTree[][] oldtrees;
	private AlloppMulLabTree mullabtree;
		
	// grjtodo this is public, copying JH - is that necessary? TreeNodeSlide accesses it.
	// 2011-06-30 parser accesses it too. 
	// Parameter or Parameter.Default ?? (a Java thing I don't get)
    public final Parameter popvalues;

	public final static boolean DBUGTUNE = false;


	public enum LegType {
		NONE, TWOBRANCH, ONEBRANCH, JOINED, NODIPLOIDS
	};

	public final static int DITREES = 0;
	public final static int TETRATREES = 1;
	public final static int NUMBEROFPLOIDYLEVELS = 2;
	
	
	
	/*
	 * Constructors. 
	 * 
	 */ 
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apspecies, double popvalue, boolean onehyb) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		apsp = apspecies;
		addModel(apsp);
		
		int ndips = apsp.SpeciesWithinPloidyLevel(2).length;
		if (onehyb) {
			if (ndips == 0) {
				makeInitialOneTetraTreeNetwork();	
			} else if (ndips == 2) {
				makeInitialOneTetraTreeTwoDiploidsNetwork(LegType.JOINED); // soon ood 2012-05-07
			} else if (ndips > 2){
				makeInitialNTetsNDipsNetwork();
			} else {
				assert false; // ndips == 1 not allowed
			}
		} else {
			assert false; // grjtodo morethanonetree need to deal with other cases
		}
		

		double maxrootheight = 0.0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			double height = trees[DITREES][i].getMaxHeight();
			if (height > maxrootheight) { maxrootheight = height; }
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			double height = trees[TETRATREES][i].getMaxHeight();
			if (height > maxrootheight) { maxrootheight = height; }
		}
		double scale = 0.99 * apsp.initialMinGeneNodeHeight() / maxrootheight;
		scaleAllHeights(scale);
		
		int npopparams = numberOfPopParameters();
		this.popvalues = new Parameter.Default(npopparams, popvalue);
		addVariable(popvalues);
		
		mullabtree = new AlloppMulLabTree(trees, apsp, popvalues);
		
        Logger.getLogger("dr.evomodel.speciation.allopolyploid").info("\tConstructing an allopolyploid network,  please cite:\n"
                + Citable.Utils.getCitationString(this));
	}	
	
		
	
	/*
	 * This constructor is for testing.
	 */
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apsp, 
			AlloppSpeciesNetworkModelTEST.NetworkConversionTEST nmltTEST) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		this.apsp = apsp;
		popvalues = null;
	}
	
	
	/*
	 * This constructor is for testing.
	 */
	public AlloppSpeciesNetworkModel(AlloppSpeciesBindings testASB,
			AlloppSpeciesNetworkModelTEST.LogLhoodGTreeInNetworkTEST llgtnTEST) {
		super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
		apsp = testASB;
		
		makeInitialOneTetraTreeNetwork(llgtnTEST);			
		
		this.popvalues = new Parameter.Default(llgtnTEST.popvalues);
		
		mullabtree = new AlloppMulLabTree(trees, apsp, popvalues);
	}

	
	

	
	
	public List<Citation> getCitations() {
		List<Citation> citations = new ArrayList<Citation>();
		citations.add(new Citation(
				new Author[]{
						new Author("GR", "Jones")
				},
				Citation.Status.IN_PREPARATION
		));
		return citations;
	}
	
	
	public String toString() {
		int ngt = apsp.numberOfGeneTrees();
		String nl = System.getProperty("line.separator");
		String s = nl + mullabtree.asText() + nl;
		for (int g = 0; g < ngt; g++) {
			s += "Gene tree " + g + nl;
			s += apsp.genetreeAsText(g) + nl;
			s += apsp.seqassignsAsText(g) + nl;
		}
		s += nl;
		return s;
	}
	
	
	
	public LogColumn[] getColumns() {
		LogColumn[] columns = new LogColumn[1];
		columns[0] = new LogColumn.Default("    MUL-tree and gene trees", this);
		return columns;
	}

	
	/*
	 * Stretches or squashes the whole network. Used by constructors and
	 * MCMC operators.
	 */
	public int scaleAllHeights(double scale) {
		int count = 0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			count += trees[DITREES][i].scaleAllHeights(scale);
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			count += trees[TETRATREES][i].scaleAllHeights(scale);
		}
		return count;
	}

	
	
	public void replaceDiploidHistory(AlloppDiploidHistory adh) {
		assert adh.diphistOK();
		for (int i = 0; i < trees[TETRATREES].length; ++i) {
			AlloppDiploidHistory.HybHistory hh = adh.extractHybHistory(i);
			trees[TETRATREES][i].replaceLegs(hh);
		}
		assert adh.diphistOK();
		trees[DITREES][0] = new AlloppLeggedTree(adh.ditreeFromDipHist());
	}
	

	/*
	 * Called from AlloppSpeciesBindings to check if a node in a gene tree
	 * is compatible with the network. 
	 */
	public boolean coalescenceIsCompatible(double height, FixedBitSet union) {
		return mullabtree.coalescenceIsCompatible(height, union);
	}
	
	
    /*
     * Called from AlloppSpeciesBindings to remove coalescent information
     * from branches of mullabtree. Required before call to recordCoalescence
     */
	public void clearCoalescences() {
		mullabtree.clearCoalescences();
	}	

	
	/*
	 * Called from AlloppSpeciesBindings to add a node from a gene tree
	 * to its branch in mullabtree.
	 */
	public void recordCoalescence(double height, FixedBitSet union) {
		mullabtree.recordCoalescence(height, union);
	}
	
	
	public void sortCoalescences() {
		mullabtree.sortCoalescences();
	}
	

	
	/*
	 * Records the number of gene lineages at nodes of mullabtree.
	 */
	public void recordLineageCounts() {
		mullabtree.recordLineageCounts();
	}	
	
	
	/*
	 * Calculates the log-likelihood for a single gene tree in the network
	 * 
	 * Requires that clearCoalescences(), recordCoalescence(), recordLineageCounts()
	 * called to fill mullabtree with information about gene tree coalescences first.
	 */
	public double geneTreeInNetworkLogLikelihood() {
		boolean noDiploids = trees[DITREES].length == 0;
		return mullabtree.geneTreeInMULTreeLogLikelihood(noDiploids);
	}	

	
	public int getTipCount() {
		int n = 0;
		for (int i = 0; i < trees[DITREES].length; i++) {
			n += trees[DITREES][i].getExternalNodeCount();
		}
		for (int i = 0; i < trees[TETRATREES].length; i++) {
			n += trees[TETRATREES][i].getExternalNodeCount();
		}
		return n;
	}
		

	public int getNumberOfDiTrees() {
		return trees[DITREES].length;
	}
	
	public int getNumberOfTetraTrees() {
		return trees[TETRATREES].length;
	}

	public int getNumberOfNodeHeightsInTree(int pl, int n) {
		return trees[pl][n].getInternalNodeCount();
	}	
	
	
	public AlloppLeggedTree getHomoploidTree(int pl, int t) {
		return trees[pl][t];
	}
	
	AlloppLeggedTree[][] getAllHomoploidTrees() {
		return trees;
	}
	
	
	 AlloppSpeciesBindings getSpeciesBindings() {
		 return apsp;
	 }
	
	
	public double getMaxFootHeight() {
		// grjtodo tetraonly
		double maxfh = 0.0;
		for (AlloppLeggedTree ttree : trees[TETRATREES]) {
			for (int i = 0;  i < ttree.getNumberOfLegs();  ++i) {
				maxfh = Math.max(maxfh, ttree.getFootHeight(i));
			}
		}
		return maxfh;
	}
	
	
	public String mullabTreeAsText() {
		return mullabtree.asText();
	}

	

	// grjtodo beginNetworkEdit(), endNetworkEdit():
	// I am doing an analogous thing to what speciesTreeModel methods
	// beginTreeEdit(), endTreeEdit() do.
	// But I don't understand the purpose.
	// 2011-07-06, OK I am starting to: remake mullabtree() after edits.
	// Could do oldmullabtree = mullabtree, etc, instead.
	// Could use dirty flags instead, and remake on demand.
    public boolean beginNetworkEdit() {
    	boolean beingEdited = false;
    	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			beingEdited = beingEdited || trees[i][j].beginTreeEdit();
    		}
    	} 
    	return beingEdited;
    }


    public void endNetworkEdit() {
     	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			trees[i][j].endTreeEdit();
    		}
    	}
     	mullabtree = new AlloppMulLabTree(trees, apsp, popvalues);
     	fireModelChanged();
    }

	

	
	public String getName() {
		return getModelName();
	}

	
    //  based on SpeciesTreeModel
    //  grjtodo internalTreeOP remaining to do: scaling without enforcing consitency
 	public int scale(double scaleFactor, int nDims) throws OperatorFailedException {
    	assert scaleFactor > 0;
    	if (nDims <= 0) {
    		beginNetworkEdit();
    		int count = scaleAllHeights(scaleFactor);
    		endNetworkEdit();
    		fireModelChanged(this, 1);
    		return count;
    	} else {
    		if (nDims != 1) {
    			throw new OperatorFailedException("not implemented for count != 1");
    		}
    		/*
            if (internalTreeOP == null) {
                internalTreeOP = new TreeNodeSlide(this, species, 1);
            }

            internalTreeOP.operateOneNode(scaleFactor);*/
    		fireModelChanged(this, 1);
    		return nDims;
    	}	
    }


    // MCMC operator which moves the legs of a tetraploid subtree, 
    // ie, changes the way it joins the diploid tree
	public void moveLegs() {
		
		
		/*
		2012-05-07 hope ood
		
		
		
		// grjtodo tetraonly morethanonetree
		// 2011-08-31 For now, the new legs are chosen from the same distribution regardless of current state.
		// I am not clear what this distribution should be. There are always two times,
		// both between the hybridizaton time and the diploid root. For the two-diploids case,
		// they are straightforward. Then there is the topology. There are five distinguishable
		// topologies in the two-diploids case: both legs left, both legs right, joined then left, 
		// joined then right, and one leg to each. I think the latter should be regarded as two cases
		// (most recent to left vs most recent to right) and it is here.
		assert trees[DITREES].length == 1;
		assert trees[DITREES][0].getExternalNodeCount() == 2;
		assert trees[TETRATREES].length == 1;
		int choice = MathUtils.nextInt(3);
		double dirooth = trees[DITREES][0].getRootHeight();
		if (choice == 0) {
			// change times but not topology or leg order
			if (MathUtils.nextInt(2) == 0) {
				trees[TETRATREES][0].moveMostRecentLegHeight();
			} else {
				trees[TETRATREES][0].moveMostAncientLegHeight(dirooth);
			}
		} else if (choice == 1) {
			// change times, maybe leg order but not topology
			trees[TETRATREES][0].moveSplitOrLeg(dirooth);
		} else {
			// change topology but not times 
			trees[TETRATREES][0].moveLegTopology(diploidtipbitset(0), diploidtipbitset(1));
		}
		
		*/
		
	}    

    
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.handleModelChangedEvent " + model.getId());
		fireModelChanged();
	}

	protected final void handleVariableChangedEvent(Variable variable,
							int index, Parameter.ChangeType type) {
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.handleVariableChangedEvent" + variable.getId());

	}

	
	
	
	protected void storeState() {
		oldtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		for (int i=0; i<trees.length; i++) {
			oldtrees[i] = new AlloppLeggedTree[trees[i].length];
		}
    	for (int i=0; i<trees.length; i++) {
    		for (int j=0; j<trees[i].length; j++) {
    			oldtrees[i][j] = new AlloppLeggedTree(trees[i][j]);
    		}
    	}
    	// addVariable(popvalues) deals with popvalues
    	
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.storeState()");
	}

	
	
	protected void restoreState() {
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		for (int i=0; i<oldtrees.length; i++) {
			trees[i] = new AlloppLeggedTree[oldtrees[i].length];
		}
    	for (int i=0; i<oldtrees.length; i++) {
    		for (int j=0; j<oldtrees[i].length; j++) {
    			trees[i][j] = new AlloppLeggedTree(oldtrees[i][j]);
    		}
    	}		
    	mullabtree = new AlloppMulLabTree(trees, apsp, popvalues);
    	// addVariable(popvalues) deals with popvalues
    	
		if (DBUGTUNE)
			System.err.println("AlloppSpeciesNetworkModel.restoreState()");
	}

	
	protected void acceptState() {
	}
	
	
	
	
	
    /*
     * For simple case of one tetraploid tree and no diploids.
     * Assume a history before root of a diploid speciating
     * at time s, the two diploids (or two descendants) forming
     * a hybrid at time h, which speciates at time r, the root
     * of the tetraploid tree.
     * 
     * heights is used for testing in cases of 2 or 3 tetraploids
     * (2011-06-13). Here heights[] supply r,h,s or d,r,h,s where
     * r,h,s are as above and d in the 3 tetraploid case is the
     * most recent tetraploid divergence.
     * 
     * The diploids go extinct or are not sampled.
     * 
     * There will be one population parameter at each tip,
     * one at the bottom of each branch within tetraploid tree
     * including the root node which goes down to h, and one
     * more population parameter for both diploids. The root
     * of the diploids will have population derived from this.
     */
	private void makeInitialOneTetraTreeNetwork() {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		trees[DITREES] = new AlloppLeggedTree[0];
		
		Taxon[] spp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(
				                       spp, LegType.NODIPLOIDS);
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
		
	}
	
	
	
	
	
    /*
     * For case of one tetraploid tree and two diploids.
     * Assume that a single diploid splits at the root,
     * and that a single hybridization event takes place between
     * the two initial diploids or descendants of them. Both initial
     * diploids leave exactly one descendant in the sample.
     * 
     * In other words, the data consist of two diploids and one
     * or more allotetraploids, and it is assumed that all
     * species arose from the MRCA of the two diploids, and that all
     * allotetraploids arose from a single hybridization event.
     * 
     * For a single tetraploid, these (plus reflections) are the 
     * possible evolutionary histories.
     * 
     * \     ||     /    \     ||     /   \     ||     /
     *  \....||..../      \....||..| /     \  |.||..| /
     *   \        /        \       |/       \ |     |/
     *    \      /          \      /         \|     /
     *     \    /            \    /           \    /
     *      \  /              \  /             \  /
     *       \/                \/               \/            
     *      
     * \  ||        /   \    ||      /
     *  \.||.|     /     \ |.||.|   /
     *   \   |    /       \|    |  /
     *    \  |   /         \   /  /
     *     \ |  /           \ |  /
     *      \| /             \| /
     *       \/               \/      
     *      
     * \    ||      /
     *  \ |.||.|   /
     *   \ \  /   /
     *    \ \/   /
     *     \/   /
     *      \  /
     *       \/           
     */
	private void makeInitialOneTetraTreeTwoDiploidsNetwork(LegType legtype) {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		Taxon[] dipspp = apsp.SpeciesWithinPloidyLevel(2);		
		AlloppLeggedTree diptree = new AlloppLeggedTree(dipspp, LegType.NONE);
		trees[DITREES] = new AlloppLeggedTree[1];
		trees[DITREES][0] = diptree;
		
		Taxon[] tetspp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(tetspp, legtype);
		// make the diploid root older than the earliest foot
		diptree.setNodeHeight(diptree.getRoot(), 1.1*tettree.getMaxHeight());
		FixedBitSet dip;
		switch (legtype) {
		case TWOBRANCH:
			// attach feet, one to each diploid
			for (int i = 0; i < 2; i++) {
				dip = diploidtipbitset(i);
				tettree.setFootUnion(i, dip);
			}
			break;
		case ONEBRANCH:
			// attach feet, both to one diploid
			dip = diploidtipbitset(0);
			for (int i = 0; i < 2; i++) {
				tettree.setFootUnion(i, dip);
			}
			break;
		case JOINED:
			// attach foot to a diploid
			dip = diploidtipbitset(0);
			tettree.setFootUnion(0, dip);
		}
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
	}
	
	
	
	
	/*
	 * Make a random initial starting network. 
	 */
	private void makeInitialNTetsNDipsNetwork() {
		// grjtodo. This is crude. Only one tet tree, never LegType.JOINED,
		// only ever joins to terminal branches, squashes tet tree to fit.
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		Taxon[] dipspp = apsp.SpeciesWithinPloidyLevel(2);
		Taxon[] tetspp = apsp.SpeciesWithinPloidyLevel(4);
		
		// make diploid tree
		AlloppLeggedTree diptree = new AlloppLeggedTree(dipspp, LegType.NONE);
		trees[DITREES] = new AlloppLeggedTree[1];
		trees[DITREES][0] = diptree;
		
		// choose random tips, decide legtype, find max height for tet tree
		int d0 = MathUtils.nextInt(dipspp.length);
		int d1 = MathUtils.nextInt(dipspp.length);
		LegType legtype = LegType.TWOBRANCH;
		if (d0 == d1) {
			legtype = LegType.ONEBRANCH;
		} else {
			legtype = LegType.TWOBRANCH;
		}
		NodeRef d0anc = diptree.getParent(diptree.getNode(d0));
		NodeRef d1anc = diptree.getParent(diptree.getNode(d1));
		double d0hgt = diptree.getNodeHeight(d0anc);	
		double d1hgt = diptree.getNodeHeight(d1anc);
		double maxteth = Math.min(d0hgt, d1hgt);
		
		// make tet tree, squash to fit into diploid tree, join up legs
		AlloppLeggedTree tettree = new AlloppLeggedTree(tetspp, legtype);
		double teth = tettree.getMaxHeight();
		tettree.scaleAllHeights(0.9 * maxteth / teth);
		FixedBitSet dip0 = diploidtipbitset(d0);
		FixedBitSet dip1 = diploidtipbitset(d1);
		tettree.setFootUnion(0, dip0);
		tettree.setFootUnion(1, dip1);
		
		// add tet tree
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;		
	}	
	
	
	
	private FixedBitSet diploidtipbitset(int i) {
		FixedBitSet dip = new FixedBitSet(apsp.numberOfSpSeqs());
		// note that ditree has been  constructed randomly, so getExternalNode(i) 
		// chooses an arbitrary node.
		String dipname = trees[DITREES][0].getNodeTaxon(trees[DITREES][0].getExternalNode(i)).getId();
		int sp = apsp.apspeciesId2index(dipname);
		int spseq = apsp.spandseq2spseqindex(sp, 0);
		// grjtodo tetraonly Don't like this way of finding the spseq index. 
		// I do similar code when constructing mullab tree, but that's OK
		// because seq is passed in. Here I am relying on diploids only 
		// having seq==0. It won't work for hexaploids joining tetras.
		// Maybe I want a apsp.firstspseqindexofspecies(species name)
		dip.set(spseq);
		return dip;
	}
	

	
	/*
	 * 2011-08-04. For no diploids, I am using one population parameter for
	 * of diploid history (1), one at tips of tetras (n), one at all branch bottoms
	 * in tettree (2n-2), and one for post-hybridization (1). Total 3n
	 * 
	 * For d>1 diploids, 3d-2 for ditree, and 3t+1 for a tetratree with t tips.
	 * 
	 * ditree is:  d (tips) + 2d-2 ( branch bottoms) = 3d-2
	 * 
	 * A tetratree is: t (tips) + 2t-2 ( branch bottoms) + 
	 *             2 (feet or foot+split) + 1 (hybridization) = 3t+1
	 * 
	 * Not sure what I want for other scenarios. 
	 */
	private int numberOfPopParameters() {
		int dim = 0;
		if (trees[DITREES].length == 0) {
			assert trees[TETRATREES].length == 1;
			int ntetratips = trees[TETRATREES][0].getExternalNodeCount();
			dim = 3*ntetratips; 
		} else {
			int nditips = trees[DITREES][0].getExternalNodeCount();
			if (nditips > 1) {
				dim += 3*nditips - 2;
				for (int i = 0; i < trees[TETRATREES].length; i++) {
					int ntetratips = trees[TETRATREES][i].getExternalNodeCount();
					dim += 3*ntetratips + 1;
				}
			} else {
				assert false; // grjtodo useful case??
			}
		}
		return dim;
	}


	

	public int getTaxonCount() {
		return mullabtree.simptree.getTaxonCount();
	}




	public Taxon getTaxon(int taxonIndex) {
		return mullabtree.simptree.getTaxon(taxonIndex);
	}




	public String getTaxonId(int taxonIndex) {
		return mullabtree.simptree.getTaxonId(taxonIndex);
	}




	public int getTaxonIndex(String id) {
		return mullabtree.simptree.getTaxonIndex(id);
	}




	public int getTaxonIndex(Taxon taxon) {
		// grjtodo. i presume references are unique
		return mullabtree.simptree.getTaxonIndex(taxon);
	}




	public List<Taxon> asList() {
		return mullabtree.simptree.asList();
	}




	public Object getTaxonAttribute(int taxonIndex, String name) {
		return mullabtree.simptree.getTaxonAttribute(taxonIndex, name);
	}




	public String getId() {
		return mullabtree.simptree.getId();
	}




	public void setId(String id) {
		mullabtree.simptree.setId(id);
		
	}




	public Iterator<Taxon> iterator() {
		return mullabtree.simptree.iterator();
	}




	public Type getUnits() {
		return mullabtree.simptree.getUnits();
	}




	public void setUnits(Type units) {
		mullabtree.simptree.setUnits(units);
		
	}




	public void setAttribute(String name, Object value) {
		mullabtree.simptree.setAttribute(name, value);
		
	}




	public Object getAttribute(String name) {
		return mullabtree.simptree.getAttribute(name);
	}




	public Iterator<String> getAttributeNames() {
		return mullabtree.simptree.getAttributeNames();
	}




	public NodeRef getRoot() {
		return mullabtree.simptree.getRoot();
	}




	public int getNodeCount() {
		return mullabtree.simptree.getNodeCount();
	}




	public NodeRef getNode(int i) {
		return mullabtree.simptree.getNode(i);
	}




	public NodeRef getInternalNode(int i) {
		return mullabtree.simptree.getInternalNode(i);
	}




	public NodeRef getExternalNode(int i) {
		return mullabtree.simptree.getExternalNode(i);
	}




	public int getExternalNodeCount() {
		return mullabtree.simptree.getExternalNodeCount();
	}




	public int getInternalNodeCount() {
		return mullabtree.simptree.getInternalNodeCount();
	}




	public Taxon getNodeTaxon(NodeRef node) {
		return mullabtree.simptree.getNodeTaxon(node);
	}




	public boolean hasNodeHeights() {
		return true;
	}




	public double getNodeHeight(NodeRef node) {
		return mullabtree.simptree.getNodeHeight(node);
	}




	public boolean hasBranchLengths() {
		return true;
	}




	public double getBranchLength(NodeRef node) {
		return mullabtree.simptree.getBranchLength(node);
	}




	public double getNodeRate(NodeRef node) {
		return mullabtree.simptree.getNodeRate(node);
	}




	public Object getNodeAttribute(NodeRef node, String name) {
		return mullabtree.simptree.getNodeAttribute(node, name);
	}




	public Iterator getNodeAttributeNames(NodeRef node) {
		return mullabtree.simptree.getNodeAttributeNames(node);
	}




	public boolean isExternal(NodeRef node) {
		return mullabtree.simptree.isExternal(node);
	}




	public boolean isRoot(NodeRef node) {
		return mullabtree.simptree.isRoot(node);
	}




	public int getChildCount(NodeRef node) {
		int cc = mullabtree.simptree.getChildCount(node);
		assert cc == 2;
		return cc;
	}




	public NodeRef getChild(NodeRef node, int j) {
		return mullabtree.simptree.getChild(node, j);
	}




	public NodeRef getParent(NodeRef node) {
		return mullabtree.simptree.getParent(node);
	}




	public Tree getCopy() {
		return mullabtree.simptree.getCopy();
	}


	public boolean logNow(long state) {
		return true;
	}    		

	

	
	
/* *********************** TEST CODE **********************************/	
	
	
	

	
	/*
	 * Test of conversion from network to mullab tree
	 * 	 * 2011-05-07 It is called from testAlloppSpeciesNetworkModel.java.
	 * I don't know how to put the code in there without
	 * making lots public here.
	 */
	// grjtodo. should be possible to pass stuff in nmltTEST. Currently
	// it just signals that this is indeed a test.

	

	public String testExampleNetworkToMulLabTree(
			AlloppSpeciesNetworkModelTEST.NetworkConversionTEST nmltTEST) {
		AlloppLeggedTree[][] testtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		int ntaxa = apsp.numberOfSpecies();
		Taxon[] spp = new Taxon[ntaxa];
		for (int tx = 0; tx < ntaxa; ++tx) {
			spp[tx] = new Taxon(apsp.apspeciesName(tx));
		}
		// 1,2,3 (names b,c,d) are tets, 0,4 are dips (names a,e)

		double tetheight = 0.0;
		// case 1. one tettree with one foot in each diploid branch
		// case 2. one tettree with both feet in one diploid branch
		// case 3. one tettree with one joined
		// case 4. two tettrees, 2+1, first with one foot in each diploid
		// branch, second joined
		// case 5. three tettrees, 1+1+1, one of each type of feet, as in cases 1-3

		int ntettrees = 0;
		switch (nmltTEST.testcase) {
		case 1:
		case 2:
		case 3:
			ntettrees = 1;
			break;
		case 4:
			ntettrees = 2;
			break;
		case 5:
			ntettrees = 3;
			break;
		}
		AlloppLeggedTree tettrees[] = new AlloppLeggedTree[ntettrees];

		Taxon[] tets123 = { spp[1], spp[2], spp[3] };
		Taxon[] tets12 = { spp[1], spp[2] };
		Taxon[] tets1 = { spp[1] };
		Taxon[] tets2 = { spp[2] };
		Taxon[] tets3 = { spp[3] };
		switch (nmltTEST.testcase) {
		case 1:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.TWOBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 2:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.ONEBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 3:
			tettrees[0] = new AlloppLeggedTree(tets123, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 4:
			tettrees[0] = new AlloppLeggedTree(tets12, nmltTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets3, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			break;
		case 5:
			tettrees[0] = new AlloppLeggedTree(tets1, nmltTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets2, nmltTEST, LegType.ONEBRANCH, 0.0);
			tettrees[2] = new AlloppLeggedTree(tets3, nmltTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			tetheight = Math.max(tetheight, tettrees[2].getMaxFootHeight());
			break;
		}

		AlloppLeggedTree ditrees[] = new AlloppLeggedTree[1];
		Taxon[] dips = { spp[0], spp[4] };
		ditrees[0] = new AlloppLeggedTree(dips, nmltTEST, LegType.NONE, tetheight + 4.0);

		testtrees[0] = ditrees;
		testtrees[1] = tettrees;


		FixedBitSet a = new FixedBitSet(8);
		a.set(0);
		FixedBitSet e = new FixedBitSet(8);
		e.set(7);
		switch (nmltTEST.testcase) {
		case 1:
			// leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			break;
		case 2:
			// both legs to node a,
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, a);
			break;
		case 3:
			// only leg to node a
			testtrees[TETRATREES][0].setFootUnion(0, a);
			break;
		case 4:
			// first tet tree (with two tips): leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			// second tet tree, only leg to node a
			testtrees[TETRATREES][1].setFootUnion(0, a);
			break;
		case 5:
			// first tet tree. leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, e);
			// second tet tree. both legs to node a,
			testtrees[TETRATREES][1].setFootUnion(0, a);
			testtrees[TETRATREES][1].setFootUnion(1, a);
			// third tet tree. only leg to node a
			testtrees[TETRATREES][2].setFootUnion(0, a);
			break;
		}
		AlloppMulLabTree testmullabtree = new AlloppMulLabTree(testtrees, apsp, popvalues);
		String newick = testmullabtree.mullabTreeAsNewick();
		return newick;
	}

	
	
	/*
	 * Test of conversion from network to diploid history
	 */
	// grjtodo implementation issues same as testExampleNetworkToDipHist()	
	public AlloppDiploidHistory testExampleNetworkToDipHist2d3t(
			AlloppSpeciesNetworkModelTEST.NetworkConversionTEST netconvTEST) {
		AlloppLeggedTree[][] testtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		int ntaxa = apsp.numberOfSpecies();
		Taxon[] spp = new Taxon[ntaxa];
		for (int tx = 0; tx < ntaxa; ++tx) {
			spp[tx] = new Taxon(apsp.apspeciesName(tx));
		}
		// 1,2,3 (names z,y,x) are tets, 0,4 are dips (names a,b)

		double tetheight = 0.0;
		// case 1. one tettree with one foot in each diploid branch
		// case 2. one tettree with both feet in one diploid branch
		// case 3. one tettree with one joined
		// case 4. two tettrees, 2+1, first with one foot in each diploid
		// branch, second joined
		// case 5. three tettrees, 1+1+1, one of each type of feet, as in cases 1-3

		int ntettrees = 0;
		switch (netconvTEST.testcase) {
		case 1:
		case 2:
		case 3:
			ntettrees = 1;
			break;
		case 4:
			ntettrees = 2;
			break;
		case 5:
			ntettrees = 3;
			break;
		}
		AlloppLeggedTree tettrees[] = new AlloppLeggedTree[ntettrees];

		Taxon[] tets123 = { spp[1], spp[2], spp[3] };
		Taxon[] tets12 = { spp[1], spp[2] };
		Taxon[] tets1 = { spp[1] };
		Taxon[] tets2 = { spp[2] };
		Taxon[] tets3 = { spp[3] };
		switch (netconvTEST.testcase) {
		case 1:
			tettrees[0] = new AlloppLeggedTree(tets123, netconvTEST, LegType.TWOBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 2:
			tettrees[0] = new AlloppLeggedTree(tets123, netconvTEST, LegType.ONEBRANCH, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 3:
			tettrees[0] = new AlloppLeggedTree(tets123, netconvTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			break;
		case 4:
			tettrees[0] = new AlloppLeggedTree(tets12, netconvTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets3, netconvTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			break;
		case 5:
			tettrees[0] = new AlloppLeggedTree(tets1, netconvTEST, LegType.TWOBRANCH, 0.0);
			tettrees[1] = new AlloppLeggedTree(tets2, netconvTEST, LegType.ONEBRANCH, 0.0);
			tettrees[2] = new AlloppLeggedTree(tets3, netconvTEST, LegType.JOINED, 0.0);
			tetheight = tettrees[0].getMaxFootHeight();
			tetheight = Math.max(tetheight, tettrees[1].getMaxFootHeight());
			tetheight = Math.max(tetheight, tettrees[2].getMaxFootHeight());
			break;
		}

		AlloppLeggedTree ditrees[] = new AlloppLeggedTree[1];
		Taxon[] dips = { spp[0], spp[4] };
		ditrees[0] = new AlloppLeggedTree(dips, netconvTEST, LegType.NONE, tetheight + 4.0);

		testtrees[0] = ditrees;
		testtrees[1] = tettrees;


		FixedBitSet a = new FixedBitSet(8);
		a.set(0);
		FixedBitSet b = new FixedBitSet(8);
		b.set(7);
		switch (netconvTEST.testcase) {
		case 1:
			// leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, b);
			break;
		case 2:
			// both legs to node a,
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, a);
			break;
		case 3:
			// only leg to node a
			testtrees[TETRATREES][0].setFootUnion(0, a);
			break;
		case 4:
			// first tet tree (with two tips): leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, b);
			// second tet tree, only leg to node a
			testtrees[TETRATREES][1].setFootUnion(0, a);
			break;
		case 5:
			// first tet tree. leg 0 to node a, leg to node e
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, b);
			// second tet tree. both legs to node a,
			testtrees[TETRATREES][1].setFootUnion(0, a);
			testtrees[TETRATREES][1].setFootUnion(1, a);
			// third tet tree. only leg to node a
			testtrees[TETRATREES][2].setFootUnion(0, a);
			break;
		}
		AlloppDiploidHistory testdiphist = new AlloppDiploidHistory(testtrees, apsp);
		return testdiphist;
	}

	
	
	/*
	 * Test of conversion from network to diploid history
	 */
	// grjtodo implementation issues same as testExampleNetworkToDipHist()	
	public AlloppDiploidHistory testExampleNetworkToDipHist3d1t(
			AlloppSpeciesNetworkModelTEST.NetworkConversionTEST netconvTEST) {
		AlloppLeggedTree[][] testtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		int ntaxa = apsp.numberOfSpecies();
		Taxon[] spp = new Taxon[ntaxa];
		for (int tx = 0; tx < ntaxa; ++tx) {
			spp[tx] = new Taxon(apsp.apspeciesName(tx));
		}
		// 0,1,2 (names a,b,c) are dips, 3 is tet (name z)

		FixedBitSet a = new FixedBitSet(5);
		a.set(0);
		FixedBitSet b = new FixedBitSet(5);
		b.set(1);
		FixedBitSet c = new FixedBitSet(5);
		c.set(2);
		FixedBitSet ab = new FixedBitSet(5);
		ab.set(0);
		ab.set(1);
		
		// case 1. ((a,b),c), tet from ab and c
		// case 2. ((a,b),c), tet from a and a
		// case 3. ((a,b),c), tet from ab split

		int ntettrees = 1;
		AlloppLeggedTree tettrees[] = new AlloppLeggedTree[ntettrees];

		Taxon[] tets = { spp[3] };
		switch (netconvTEST.testcase) {
		case 1:
			tettrees[0] = new AlloppLeggedTree(tets, netconvTEST, LegType.TWOBRANCH, 0.0);
			break;
		case 2:
			tettrees[0] = new AlloppLeggedTree(tets, netconvTEST, LegType.ONEBRANCH, 0.0);
			break;
		case 3:
			tettrees[0] = new AlloppLeggedTree(tets, netconvTEST, LegType.JOINED, 0.0);
			break;
		}

		AlloppLeggedTree ditrees[] = new AlloppLeggedTree[1];
		Taxon[] dips = { spp[0], spp[1], spp[2] };
		ditrees[0] = new AlloppLeggedTree(dips, netconvTEST, LegType.NONE, 0.0);

		testtrees[0] = ditrees;
		testtrees[1] = tettrees;

		switch (netconvTEST.testcase) {
		case 1:
			// leg 0 to node ab, leg 1 to node c
			testtrees[TETRATREES][0].setFootUnion(0, ab);
			testtrees[TETRATREES][0].setFootUnion(1, c);
			break;
		case 2:
			// legs 0 and 1 to node a
			testtrees[TETRATREES][0].setFootUnion(0, a);
			testtrees[TETRATREES][0].setFootUnion(1, a);
			break;
		case 3:
			// leg 0 to node ab, no leg 1
			testtrees[TETRATREES][0].setFootUnion(0, ab);
			break;
		}
		AlloppDiploidHistory testdiphist = new AlloppDiploidHistory(testtrees, apsp);
		return testdiphist;
	}
	
	
	
	
	/*
	 * Test of conversion from network to diploid history
	 */
	// grjtodo implementation issues same as testExampleNetworkToDipHist()	
	public AlloppDiploidHistory testExampleNetworkToDipHist4d2t(
			AlloppSpeciesNetworkModelTEST.NetworkConversionTEST netconvTEST) {
		AlloppLeggedTree[][] testtrees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		
		int ntaxa = apsp.numberOfSpecies();
		Taxon[] spp = new Taxon[ntaxa];
		for (int tx = 0; tx < ntaxa; ++tx) {
			spp[tx] = new Taxon(apsp.apspeciesName(tx));
		}
		assert(ntaxa==6);
		// 0,1,2,3 (names a,b,c,d) are dips, 4,5 are tets (names z,y)
		
		FixedBitSet a =   new FixedBitSet(8);   a.set(0);
		FixedBitSet b =   new FixedBitSet(8);   b.set(1);
		FixedBitSet c =   new FixedBitSet(8);   c.set(2);
		FixedBitSet d =   new FixedBitSet(8);   d.set(3);
		FixedBitSet ab =  new FixedBitSet(8);  ab.set(0);  ab.set(1);
		FixedBitSet abc = new FixedBitSet(8); abc.set(0); abc.set(1); abc.set(2);
		
		// case 1. (((a,b),c),d), tet from ab and d
		int ntettrees = 1;
		AlloppLeggedTree tettrees[] = new AlloppLeggedTree[ntettrees];

		Taxon[] tets = { spp[4], spp[5] };
		switch (netconvTEST.testcase) {
		case 1:
			tettrees[0] = new AlloppLeggedTree(tets, netconvTEST, LegType.TWOBRANCH, 0.0);
			break;
		}

		AlloppLeggedTree ditrees[] = new AlloppLeggedTree[1];
		Taxon[] dips = { spp[0], spp[1], spp[2], spp[3] };
		ditrees[0] = new AlloppLeggedTree(dips, netconvTEST, LegType.NONE, 0.0);

		testtrees[0] = ditrees;
		testtrees[1] = tettrees;

		switch (netconvTEST.testcase) {
		case 1:
			// leg 0 to node ab, leg 1 to node d
			testtrees[TETRATREES][0].setFootUnion(0, ab);
			testtrees[TETRATREES][0].setFootUnion(1, d);
			break;
		}
		AlloppDiploidHistory testdiphist = new AlloppDiploidHistory(testtrees, apsp);
		return testdiphist;
	}

	
	/*
	 * for testing
	 */
	private void makeInitialOneTetraTreeNetwork(AlloppSpeciesNetworkModelTEST.LogLhoodGTreeInNetworkTEST llgtnTEST) {
		
		trees = new AlloppLeggedTree[NUMBEROFPLOIDYLEVELS][];
		trees[DITREES] = new AlloppLeggedTree[0];

		Taxon[] spp = apsp.SpeciesWithinPloidyLevel(4);
		AlloppLeggedTree tettree = new AlloppLeggedTree(
				                       spp, LegType.NODIPLOIDS, llgtnTEST);
		trees[TETRATREES] = new AlloppLeggedTree[1];
		trees[TETRATREES][0] = tettree;
		
	}

	
	

}
