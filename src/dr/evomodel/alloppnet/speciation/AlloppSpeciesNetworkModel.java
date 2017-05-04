/*
 * AlloppSpeciesNetworkModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.alloppnet.speciation;



import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.alloppnet.parsers.AlloppSpeciesNetworkModelParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.Scalable;
import dr.math.MathUtils;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import jebl.util.FixedBitSet;
import java.util.*;
import java.util.logging.Logger;


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
//  grjtodo-oneday JH's SpeciesTreeModel implements
// MutableTree, TreeTraitProvider, TreeLogger.LogUpon, Scalable
// not clear how much of those sensible here.
// AlloppLeggedTree implements MutableTree, TreeLogger.LogUpon.
// Nothing so far does TreeTraitProvider.
public class AlloppSpeciesNetworkModel extends AbstractModel implements
        Scalable, Units, Citable, Tree, TreeTraitProvider, TreeLogger.LogUpon {

    private final AlloppSpeciesBindings apsp;
    private AlloppDiploidHistory adhist;
    private AlloppDiploidHistory oldadhist;
    private ArrayList<AlloppLeggedTree> tettrees;
    private ArrayList<AlloppLeggedTree> oldtettrees;
    private int nofdiploids;
    private int noftetraploids;
    private boolean onehybridization;
    private boolean diploidrootisroot;
    private TreeTrait tti;
    private TreeTrait hh;

    private AlloppMulLabTree mullabtree;

    private ParametricDistributionModel hybridPopModel;

    // The Parameters are public, copying JH - is that necessary? TreeNodeSlide accesses it.
    // 2011-06-30 parser accesses it too.
    // Parameter or Parameter.Default ?? (a Java thing I don't get)
    public final Parameter tippopvalues;
    public final Parameter rootpopvalues;

    public final Parameter logginghybpopvalues;
    private double [] hybpopvalues;
    private double [] oldhybpopvalues;


    public final static boolean DBUGTUNE = false;




    private class TetTreeIndexTrait implements TreeTrait<String> {
        TetTreeIndexTrait() {}

        public String getTraitName() {
            return "tti";
        }

        public TreeTrait.Intent getIntent() {
            return Intent.BRANCH;
        }

        @Override
        public Class getTraitClass() {
            return String.class;
        }

        public String getTrait(Tree tree, NodeRef node) {
            assert tree == mullabtree;
            return (String)getNodeAttribute(node, "tti");
        }

        @Override
        public String getTraitString(Tree tree, NodeRef node) {
            return "" + getNodeAttribute(node, "tti");
        }

        @Override
        public boolean getLoggable() {
            return true;
        }
    }



    // grjtodo-soon hybheights in TreeAnnotator is not working as hoped.
    // Nodes may or may not
    // have hybheights, so means, medians, etc, combine -1 with valid values
    // Omitting attributes from nodes results in some `null's in the tree log file
    // which makes TreeAnnotator treat all values as discrete (as a set).
    // For now, just don't add hh to list in getTreeTraits()
    private class HybHeightTrait implements TreeTrait<Double> {
        HybHeightTrait() {}

        public String getTraitName() {
            return "hybhgt";
        }

        public TreeTrait.Intent getIntent() {
            return TreeTrait.Intent.NODE;
        }

        @Override
        public Class getTraitClass() {
            return Double.class;
        }

        public Double getTrait(Tree tree, NodeRef node) {
            assert tree == mullabtree;
            return (Double)getNodeAttribute(node, "hybhgt");
        }

        @Override
        public String getTraitString(Tree tree, NodeRef node) {
            return "" + getNodeAttribute(node, "hybhgt");
        }

        @Override
        public boolean getLoggable() {
            return true;
        }
    }




    /*
      * Constructors.
      *
      */
    public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apspecies,
                                     double tippopvalue, double rootpopvalue, double hybpopvalue,
                                     boolean onehyb, boolean diprootisroot) {
        super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
        apsp = apspecies;
        addModel(apsp);
        tettrees = new ArrayList<AlloppLeggedTree>();
        Taxon[] dipspp = apsp.SpeciesWithinPloidyLevel(2);
        nofdiploids = dipspp.length;
        Taxon[] tetspp = apsp.SpeciesWithinPloidyLevel(4);
        noftetraploids = tetspp.length;
        onehybridization = onehyb;
        diploidrootisroot = diprootisroot;
        makeInitialNDipsNTetsNetwork(dipspp, tetspp);

        double maxrootheight = adhist.getRootHeight();
        for (int i = 0; i < tettrees.size(); i++) {
            double height = tettrees.get(i).getRootHeight();
            if (height > maxrootheight) { maxrootheight = height; }
        }
        double scale = 0.99 * apsp.initialMinGeneNodeHeight() / maxrootheight;
        scaleAllHeights(scale);

        int ntippopparams = numberOfTipPopParameters();
        int nrootpopparams = numberOfRootPopParameters();
        int maxnhybpopparams = maxNumberOfHybPopParameters();
        tippopvalues = new Parameter.Default(ntippopparams, tippopvalue);
        rootpopvalues = new Parameter.Default(nrootpopparams, rootpopvalue);
        addVariable(tippopvalues);
        addVariable(rootpopvalues);
        // hybridization pop sizes have to be done differently because they change in number.
        hybpopvalues = new double[maxnhybpopparams];
        for (int hp = 0; hp < hybpopvalues.length; hp++) {
            hybpopvalues[hp] = hybpopvalue;
        }

        logginghybpopvalues = new Parameter.Default(hybpopvalues);
        makeLoggingHybPopParam();

        mullabtree = new AlloppMulLabTree(adhist, tettrees, apsp, tippopvalues, rootpopvalues, hybpopvalues);

        tti = new TetTreeIndexTrait();
        hh = new HybHeightTrait();

        Logger.getLogger("dr.evomodel.speciation.allopolyploid").info("\tConstructing an allopolyploid network,  please cite:\n"
                + Citable.Utils.getCitationString(this));
    }



    /*
      * This (partial) constructor is for testing conversion network to multree.
      * Real work done by testExampleNetworkToMulLabTree()
      */
    public AlloppSpeciesNetworkModel(AlloppSpeciesBindings apsp) {
        super(AlloppSpeciesNetworkModelParser.ALLOPPSPECIESNETWORK);
        this.apsp = apsp;
        tippopvalues = null;
        rootpopvalues = null;
        hybpopvalues = null;
        logginghybpopvalues = null;
    }

    // This is called from AlloppNetworkPrior (which is created after network)
    // to supply a ParametricDistributionModel for the prior on the
    // hybrid population values. It completes the construction of the network.
    public void setHybPopModel(ParametricDistributionModel pdm) {
        hybridPopModel = pdm;
    }


    /***********************************************************************************/


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SPECIES_MODELS;
    }

    @Override
    public String getDescription() {
        return "Allopolyploid Species Networks";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(
                        new Author[]{
                                new Author("Graham", "Jones"),
                                new Author("Serik", "Sagitov"),
                                new Author("Bengt", "Oxelman")
                        },
                        "Statistical Inference of Allopolyploid Species Networks in the Presence of Incomplete Lineage Sorting",
                        2013,
                        "Systematic Biology",
                        62,
                        467,
                        478,
                        Citation.Status.PUBLISHED
                ));
    }


    public boolean alloppspeciesnetworkOK() {
        for (AlloppLeggedTree tettree : tettrees) {
            if (!tettree.leggedtreeOK()) {
                return false;
            }
        }
        for (int tt = 0; tt <tettrees.size(); tt++) {
            AlloppLeggedTree tettree = getTetraploidTree(tt);
            int lftleg = tettree.getDiphistLftLeg();
            int rgtleg = tettree.getDiphistRgtLeg();
            if  (AlloppDiploidHistory.LegLorR.left != adhist.getNodeLeg(lftleg)) {
                return false;
            }
            if (tt != adhist.getNodeTettree(lftleg)) {
                return false;
            }
            if (AlloppDiploidHistory.LegLorR.right != adhist.getNodeLeg(rgtleg)) {
                return false;
            }
            if  (tt != adhist.getNodeTettree(rgtleg)) {
                return false;
            }
        }
        if (!adhist.diphistOK(diploidrootisroot)) {
            return false;
        }

        if (!mullabtree.mullabtreeOK()) {
            return false;
        }
        return true;
    }



    // for testing
    String mullabTreeAsText() {
        return mullabtree.asText();
    }




    // AbstractModel implementation

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (DBUGTUNE)
            System.err.println("AlloppSpeciesNetworkModel.handleModelChangedEvent " + model.getId());
        fireModelChanged();
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable,
                                                    int index, Parameter.ChangeType type) {
        if (DBUGTUNE)
            System.err.println("AlloppSpeciesNetworkModel.handleVariableChangedEvent" + variable.getId());

    }


    @Override
    protected void storeState() {
        oldtettrees = new ArrayList<AlloppLeggedTree>();
        for (int j=0; j<tettrees.size(); j++) {
            oldtettrees.add(new AlloppLeggedTree(tettrees.get(j)));
        }
        oldadhist = new AlloppDiploidHistory(adhist);
        oldhybpopvalues = new double[hybpopvalues.length];
        for (int i = 0; i < oldhybpopvalues.length; i++) {
            oldhybpopvalues[i] = hybpopvalues[i];
        }
        // addVariable(tippopvalues), addVariable(rootpopvalues) deal with other popvalues
        if (DBUGTUNE)
            System.err.println("AlloppSpeciesNetworkModel.storeState()");
    }


    @Override
    protected void restoreState() {
        tettrees = new ArrayList<AlloppLeggedTree>();
        for (int j=0; j<oldtettrees.size(); j++) {
            tettrees.add(new AlloppLeggedTree(oldtettrees.get(j)));
        }
        adhist = new AlloppDiploidHistory(oldadhist);
        hybpopvalues = new double[oldhybpopvalues.length];
        for (int i = 0; i < hybpopvalues.length; i++) {
            hybpopvalues[i] = oldhybpopvalues[i];
        }
        makeLoggingHybPopParam();
        // addVariable(tippopvalues), addVariable(rootpopvalues) deal with other popvalues
        mullabtree = new AlloppMulLabTree(adhist, tettrees, apsp, tippopvalues, rootpopvalues, hybpopvalues);
        if (DBUGTUNE)
            System.err.println("AlloppSpeciesNetworkModel.restoreState()");
    }


    @Override
    protected void acceptState() {
    }


    @Override
    public String toString() {
        int ngt = apsp.numberOfGeneTrees();
        String nl = System.getProperty("line.separator");

        String s = nl + adhist.asText() + nl;
        s += "noftettrees " + tettrees.size() + nl;
        for (int tt = 0; tt < tettrees.size(); tt++) {
            s += tettrees.get(tt).asText(tt);
        }
        s += mullabtree.asText() + nl;
        for (int g = 0; g < ngt; g++) {
            s += apsp.genetreeAsText(g);
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



    /****** Next bunch of methods used by MCMC operators  *****************************/


    public boolean beginNetworkEdit() {
        assert alloppspeciesnetworkOK();
        boolean beingEdited = false;
        return beingEdited;
    }


    // remake mullabtree and hybpop param  after edits.
    public void endNetworkEdit() {
        makeLoggingHybPopParam();
        mullabtree = new AlloppMulLabTree(adhist, tettrees, apsp, tippopvalues, rootpopvalues, hybpopvalues);
        assert alloppspeciesnetworkOK();
        fireModelChanged();
    }


    public boolean netAndGTreesAreCompatible() {
        for (int i = 0; i < apsp.numberOfGeneTrees(); i++) {
            if (!apsp.geneTreeFitsInNetwork(i, this)) {
                return false;
            }
        }
        return true;
    }

    // Scalable implementation
    @Override
    public String getName() {
        return getModelName();
    }


    // Scalable implementation. Stretches/squeezes whole network.
    @Override
    public int scale(double scaleFactor, int nDims, boolean testBounds) {
        assert scaleFactor > 0;
        assert nDims <= 0;
        if (nDims <= 0) {
            beginNetworkEdit();
            int count = 0;
            count += scaleAllHeights(scaleFactor);
            count += scaleAllPopValues(scaleFactor);
            endNetworkEdit();
            fireModelChanged(this, 1);
            return count;
        } else {
            //  grjtodo-oneday JH also has a internalTreeOP for nDims==1 case
            if (nDims != 1) {
                throw new UnsupportedOperationException("not implemented for count != 1");
            }
            fireModelChanged(this, 1);
            return nDims;
        }
    }

    @Override
    public boolean testBounds() {
        return true;
    }

    // finds the union of a tip (diploid or hyb-tip) in the diploid history.
    // Used by move that slides node in diploid history to check gene-tree compatibility.
    public FixedBitSet calculateDipHistTipUnion(NodeRef node) {
        if (node == null) {
            System.out.println("BUG in calculateDipHistTipUnion()");
        }
        assert node != null;
        int tt = adhist.getNodeTettree(node.getNumber());
        AlloppDiploidHistory.LegLorR leg = adhist.getNodeLeg(node.getNumber());
        int seq = (leg == AlloppDiploidHistory.LegLorR.left) ? 0 : 1;
        FixedBitSet union;
        if (tt < 0) { // ordinary tip
            union = apsp.taxonseqToTipUnion(adhist.getSlidableNodeTaxon(node), 0);
        } else {
            union = unionOfWholeTetTree(tt, seq);
        }
        return union;
    }



    public FixedBitSet unionOfWholeTetTree(int tt, int leg) {
        // fill in tips
        tettrees.get(tt).fillinTipUnions(apsp, leg);
        // fill in rest. could just take union of all, but have this function ready
        AlloppNode tetroot = (AlloppNode)tettrees.get(tt).getSlidableRoot();
        tetroot.fillinUnionsInSubtree(apsp.numberOfSpSeqs());
        return tetroot.getUnion();
    }


    public double addHybPopParam() {
        assert tettrees.size() <= hybpopvalues.length;
        double newval = hybridPopModel.quantile(MathUtils.nextDouble());
        if (newval < 1E-10) {
            newval = 1E-10; // grjtodo-soon.
            // There is a problem. quantile(3.84E-7) with sensible param values returns 4.9E-324
            // GammaDistImpl(alpha=1,beta=6.5e-5)
        }
        hybpopvalues[tettrees.size()-1] = newval;
        return hybridPopModel.logPdf(newval);
    }


    public double removeHybPopParam() {
        assert tettrees.size() < hybpopvalues.length;
        double oldval = hybpopvalues[tettrees.size()];
        hybpopvalues[tettrees.size()] = 0.0;      // set unused dimension to impossible value
        return hybridPopModel.logPdf(oldval);
    }


    public int getNumberOfTetraTrees() {
        return tettrees.size();
    }

    public boolean getDiploidRootIsRoot() {
        return diploidrootisroot;
    }


    public boolean getOneHybridization() {
        return onehybridization;
    }

    public int getNumberOfInternalNodesInTetTree(int n) {
        return tettrees.get(n).getInternalNodeCount();
    }


    public int getNumberOfInternalNodesInDipHist() {
        return adhist.getInternalNodeCount();
    }


    public AlloppLeggedTree getTetraploidTree(int t) {
        return tettrees.get(t);
    }

    public AlloppDiploidHistory getDiploidHistory() {
        return adhist;
    }


    public int getNofDiploids() {
        return nofdiploids;
    }


    // for merge and split moves
    public void setTetTree(int oldtt, AlloppLeggedTree newttree) {
        tettrees.set(oldtt, newttree);
    }

    // for split move
    public int addTetTree(AlloppLeggedTree tettree) {
        tettrees.add(tettree);
        return tettrees.size() - 1;
    }

    // for merge move
    public void removeTetree(int tt) {
        tettrees.remove(tt);
    }


    // for move that flips all seqs of tet tree and its legs
    public void flipLegsOfTetraTree(int tt) {
        int oldlftleg = tettrees.get(tt).getDiphistLftLeg();
        int oldrgtleg = tettrees.get(tt).getDiphistRgtLeg();
        AlloppDiploidHistory.LegLorR lftLorR = adhist.getNodeLeg(oldlftleg);
        AlloppDiploidHistory.LegLorR rgtLorR = adhist.getNodeLeg(oldrgtleg);
        adhist.setNodeLeg(oldlftleg, rgtLorR);
        adhist.setNodeLeg(oldrgtleg, lftLorR);
        tettrees.get(tt).setDiphistLftLeg(oldrgtleg);
        tettrees.get(tt).setDiphistRgtLeg(oldlftleg);
    }

    public void moveLegs() {
        // ood
    }



    public int maxNumberOfHybPopParameters() {
        return apsp.SpeciesWithinPloidyLevel(4).length;
    }

    public void setOneHybPopValue(int i, double v) {
        hybpopvalues[i] = v;
    }

    /******************************** next bunch for lhood calculations **************************************/


    Parameter getTipPopValues() {
        return tippopvalues;
    }

    Parameter getRootPopValues() {
        return rootpopvalues;
    }

    public double getOneHybPopValue(int i) {
        return hybpopvalues[i];
    }

    /*
      * Called from AlloppSpeciesBindings to check if a node in a gene tree
      * is compatible with the network.
      */
    boolean coalescenceIsCompatible(double height, FixedBitSet union) {
        boolean ok = mullabtree.coalescenceIsCompatible(height, union);
        return ok;
    }


    /*
    * Called from AlloppSpeciesBindings to remove coalescent information
    * from branches of mullabtree. Required before call to recordCoalescence
    */
    void clearCoalescences() {
        mullabtree.clearCoalescences();
    }


    /*
      * Called from AlloppSpeciesBindings to add a node from a gene tree
      * to its branch in mullabtree.
      */
    void recordCoalescence(double height, FixedBitSet union) {
        mullabtree.recordCoalescence(height, union);
    }


    void sortCoalescences() {
        mullabtree.sortCoalescences();
    }


    /*
      * Records the number of gene lineages at nodes of mullabtree.
      */
    void recordLineageCounts() {
        mullabtree.recordLineageCounts();
    }


    /*
      * Calculates the log-likelihood for a single gene tree in the network
      *
      * Requires that clearCoalescences(), recordCoalescence(), recordLineageCounts()
      * called to fill mullabtree with information about gene tree coalescences first.
      */
    double geneTreeInNetworkLogLikelihood() {
        return mullabtree.geneTreeInMULTreeLogLikelihood();
    }



    /********************** for logging ********************/

    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[]{tti};
    }

    public TreeTrait getTreeTrait(String key) {
        if (key.equals(tti.getTraitName())) {
            return tti;
        }
        if (key.equals(hh.getTraitName())) {
            return hh;
        }
        throw new IllegalArgumentException();
    }


    /********************************************************************************/
    /***********************   private    *******************************************/
    /********************************************************************************/


    /*
	 * Make a random initial starting network.
	 */
    private void makeInitialNDipsNTetsNetwork(Taxon[] dipspp, Taxon[] tetspp) {
        //
        double rate = 1.0; // scale later
        assert tetspp.length > 0;
        assert dipspp.length > 1;
        ArrayList<TetraTaxonGroup> tetgps = new ArrayList<TetraTaxonGroup>();
        TetraTaxonGroup gp1 = new TetraTaxonGroup();
        if (onehybridization) {
            for (int t = 0;  t < tetspp.length; t++) {
                gp1.add(tetspp[t]);
            }
            tetgps.add(gp1);
        } else {
            // Chinese restuarant process to partition tetraploids
            gp1.add(tetspp[0]);
            tetgps.add(gp1);
            for (int t = 1;  t < tetspp.length; t++) {
                double [] pdf = new double[tetgps.size() + 1];
                for (int g = 0; g < tetgps.size(); g++) {
                    pdf[g] = tetgps.get(g).size();
                }
                pdf[tetgps.size()] = 1;
                int nextg = MathUtils.randomChoicePDF(pdf);
                if (nextg == tetgps.size()) {
                    TetraTaxonGroup newgp = new TetraTaxonGroup();
                    newgp.add(tetspp[t]);
                    tetgps.add(newgp);
                } else {
                    tetgps.get(nextg).add(tetspp[t]);
                }
            }
        }
        // Make trees for each group of tetraploids
        for (int g = 0; g < tetgps.size(); g++) {
            Taxon [] gpspp = new Taxon[tetgps.get(g).size()];
            for (int t = 0; t < tetgps.get(g).size(); t++) {
                gpspp[t] = tetgps.get(g).get(t);
            }
            AlloppLeggedTree tettree = new AlloppLeggedTree(gpspp, rate);
            tettrees.add(tettree);
        }
        // Make diploid history given tetraploid subtrees
        adhist = new AlloppDiploidHistory(dipspp, tettrees, diploidrootisroot, rate, apsp);
    }


    private class TetraTaxonGroup {
        ArrayList<Taxon> tettxs;
        TetraTaxonGroup() { tettxs = new ArrayList<Taxon>(); }
        public void add(Taxon tx) { tettxs.add(tx); }
        public Taxon get(int i) { return tettxs.get(i); }
        public int size() { return tettxs.size(); }
    }


    private void makeLoggingHybPopParam() {
        for (int i = 0; i < hybpopvalues.length; i++) {
            logginghybpopvalues.setParameterValueQuietly(i, hybpopvalues[i]);
        }
    }

    /*
    * Stretches or squashes all population values. Used by MCMC operators.
    */
    private int scaleAllPopValues(double scale) {
        int count = 0;
        for (int i = 0; i < tippopvalues.getDimension(); i++) {
            tippopvalues.setParameterValue(i, scale*tippopvalues.getParameterValue(i));
            count++;
        }
        for (int i = 0; i < rootpopvalues.getDimension(); i++) {
            rootpopvalues.setParameterValue(i, scale*rootpopvalues.getParameterValue(i));
            count++;
        }
        for (int i = 0; i < tettrees.size(); i++) {
            hybpopvalues[i] *= scale;
            count++;
        }
        return count;
    }


    /*
      * Stretches or squashes all node heights. Used by constructors and
      * MCMC operators.
      */
    private int scaleAllHeights(double scale) {
        int count = adhist.scaleAllHeights(scale);
        for (int i = 0; i < tettrees.size(); i++) {
            count += tettrees.get(i).scaleAllHeights(scale);
        }
        return count;
    }




    private int numberOfTipPopParameters() {
        int nditips = apsp.SpeciesWithinPloidyLevel(2).length;
        int ntettips = apsp.SpeciesWithinPloidyLevel(4).length;
        return nditips + ntettips;
    }

    private int numberOfRootPopParameters() {
        int nditips = apsp.SpeciesWithinPloidyLevel(2).length;
        int ntettips = apsp.SpeciesWithinPloidyLevel(4).length;
        return 2*(nditips + ntettips - 1);
    }



    /***************************************************************************/
    /****************** Delgations for Tree ******************************/
    /***************************************************************************/


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
        // can set logEvery=0 in XML for multree:
        //      <logTree id="multreeFileLog" logEvery="0" fileName="C:/U....
        // and get here for debugging
        if (state == 6696) {
            System.out.println("logNow("+state+")");
        }
        if (state <= 100) {
            return true;
        }
        if (state <= 10000) {
            return (state % 100) == 0;
        }
        return (state % 10000) == 0;
    }





/* *********************** TEST CODE **********************************/


    /*
      * Test of conversion from network to mullab tree
      * 	 * 2011-05-07 It is called from testAlloppSpeciesNetworkModel.java.
      * I don't know how to put the code in there without
      * making lots public here.
      */
    // grjtodo-oneday. should be possible to pass stuff in nmltTEST. Currently
    // it just signals that this is indeed a test.



    //AR - removing this as it creates a dependency to test.dr.* which is bad...

    public String testExampleNetworkToMulLabTree(int testcase) {

        int ntaxa = apsp.numberOfSpecies();
        Taxon[] spp = new Taxon[ntaxa];
        for (int tx = 0; tx < ntaxa; ++tx) {
            spp[tx] = new Taxon(apsp.apspeciesName(tx));
        }
        // 1,2,3 (names b,c,d) are tets, 0,4 are dips (names a,e)

        double tetheight0 = 0.0;
        double tetheight1 = 0.0;
        double tetheight2 = 0.0;
        // case 1. one tettree with one foot in each diploid branch
        // case 2. one tettree with both feet in one diploid branch
        // case 3. one tettree with one joined
        // case 4. two tettrees, 2+1, first with one foot in each diploid
        // branch, second joined
        // case 5. three tettrees, 1+1+1, one of each type of feet, as in cases 1-3

        int ntettrees = 0;
        switch (testcase) {
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
        tettrees = new ArrayList<AlloppLeggedTree>(ntettrees);

        Taxon l0 = new Taxon("L0");
        Taxon l1 = new Taxon("L1");
        Taxon l2 = new Taxon("L2");
        Taxon r0 = new Taxon("R0");
        Taxon r1 = new Taxon("R1");
        Taxon r2 = new Taxon("R2");

        Taxon[] tets123 = {spp[1], spp[2], spp[3]};
        Taxon[] tets12 = {spp[1], spp[2]};
        Taxon[] tets1 = {spp[1]};
        Taxon[] tets2 = {spp[2]};
        Taxon[] tets3 = {spp[3]};
        Taxon[] dips = new Taxon[0];
        switch (testcase) {
            case 1:
                tettrees.add(new AlloppLeggedTree(tets123));
                tetheight0 = tettrees.get(0).getRootHeight();
                dips = new Taxon[] {spp[0], l0, r0, spp[4]};
                break;
            case 2:
                tettrees.add(new AlloppLeggedTree(tets123));
                tetheight0 = tettrees.get(0).getRootHeight();
                dips = new Taxon[] {spp[0], l0, r0, spp[4]};
                break;
            case 3:
                tettrees.add(new AlloppLeggedTree(tets123));
                tetheight0 = tettrees.get(0).getRootHeight();
                dips = new Taxon[] {spp[0], l0, r0, spp[4]};
                break;
            case 4:
                tettrees.add(new AlloppLeggedTree(tets12));
                tettrees.add(new AlloppLeggedTree(tets3));
                tetheight0 = tettrees.get(0).getRootHeight();
                tetheight1 = tettrees.get(1).getRootHeight();
                dips = new Taxon[] {spp[0], l0, r0, l1, r1, spp[4]};
                break;
            case 5:
                tettrees.add(new AlloppLeggedTree(tets1));
                tettrees.add(new AlloppLeggedTree(tets2));
                tettrees.add(new AlloppLeggedTree(tets3));
                tetheight0 = tettrees.get(0).getRootHeight();
                tetheight1 = tettrees.get(1).getRootHeight();
                tetheight2 = tettrees.get(2).getRootHeight();
                dips = new Taxon[] {spp[0], l0, r0, l1, r1, l2, r2, spp[4]};
                break;
        }
        assert dips.length >= 2;
        int ndhnodes = 2*dips.length - 1;
        SimpleNode[] dhnodes = new SimpleNode[ndhnodes];
        for (int n = 0; n < ndhnodes; n++) {
            dhnodes[n] = new SimpleNode();
            if (n < dips.length) {
                dhnodes[n].setTaxon(dips[n]);
            } else {
                dhnodes[n].setTaxon(new Taxon(""));
            }
        }
        int dhroot = -1;
        switch (testcase) {
            case 1:
                dhnodes[1].setHeight(tetheight0 + 1.0);
                dhnodes[2].setHeight(tetheight0 + 1.0);
                addSimpleNodeChildren(dhnodes[4], dhnodes[0], dhnodes[1], 1.0);
                addSimpleNodeChildren(dhnodes[5], dhnodes[2], dhnodes[3], 1.0);
                addSimpleNodeChildren(dhnodes[6], dhnodes[4], dhnodes[5], 1.0);
                dhroot = 6;
                break;
            case 2:
                dhnodes[1].setHeight(tetheight0 + 1.0);
                dhnodes[2].setHeight(tetheight0 + 1.0);
                addSimpleNodeChildren(dhnodes[4], dhnodes[0], dhnodes[1], 1.0);
                addSimpleNodeChildren(dhnodes[5], dhnodes[2], dhnodes[4], 1.0);
                addSimpleNodeChildren(dhnodes[6], dhnodes[3], dhnodes[5], 1.0);
                dhroot = 6;
                break;
            case 3:
                dhnodes[1].setHeight(tetheight0 + 1.0);
                dhnodes[2].setHeight(tetheight0 + 1.0);
                addSimpleNodeChildren(dhnodes[4], dhnodes[1], dhnodes[2], 1.0);
                addSimpleNodeChildren(dhnodes[5], dhnodes[0], dhnodes[4], 1.0);
                addSimpleNodeChildren(dhnodes[6], dhnodes[3], dhnodes[5], 1.0);
                dhroot = 6;
                break;
            case 4:
                dhnodes[1].setHeight(tetheight0 + 1.0);
                dhnodes[2].setHeight(tetheight0 + 1.0);
                dhnodes[3].setHeight(tetheight1 + 1.0);
                dhnodes[4].setHeight(tetheight1 + 1.0);
                addSimpleNodeChildren(dhnodes[6], dhnodes[0], dhnodes[1], 1.0);
                addSimpleNodeChildren(dhnodes[7], dhnodes[3], dhnodes[4], 1.0);
                addSimpleNodeChildren(dhnodes[8], dhnodes[6], dhnodes[7], 1.0);
                addSimpleNodeChildren(dhnodes[9], dhnodes[2], dhnodes[5], 1.0);
                addSimpleNodeChildren(dhnodes[10], dhnodes[8], dhnodes[9], 1.0);
                dhroot = 10;
                break;
            case 5:
                dhnodes[1].setHeight(tetheight0 + 1.0);
                dhnodes[2].setHeight(tetheight0 + 1.0);
                dhnodes[3].setHeight(tetheight1 + 1.0);
                dhnodes[4].setHeight(tetheight1 + 1.0);
                dhnodes[5].setHeight(tetheight2 + 1.0);
                dhnodes[6].setHeight(tetheight2 + 1.0);
                addSimpleNodeChildren(dhnodes[8], dhnodes[0], dhnodes[1], 1.0);
                addSimpleNodeChildren(dhnodes[9], dhnodes[5], dhnodes[6], 1.0);
                addSimpleNodeChildren(dhnodes[10], dhnodes[2], dhnodes[7], 1.0);
                addSimpleNodeChildren(dhnodes[11], dhnodes[3], dhnodes[8], 1.0);
                addSimpleNodeChildren(dhnodes[12], dhnodes[4], dhnodes[11], 1.0);
                addSimpleNodeChildren(dhnodes[13], dhnodes[9], dhnodes[12], 1.0);
                addSimpleNodeChildren(dhnodes[14], dhnodes[10], dhnodes[13], 1.0);
                dhroot = 14;
                break;
        }
        AlloppDiploidHistory adhist = new AlloppDiploidHistory(dhnodes, dhroot, tettrees, true, apsp);
        int ntippopparams = numberOfTipPopParameters();
        int nrootpopparams = numberOfRootPopParameters();
        int maxnhybpopparams = maxNumberOfHybPopParameters();
        Parameter testtippopvalues = new Parameter.Default(ntippopparams);
        Parameter testrootpopvalues = new Parameter.Default(nrootpopparams);
        double [] testhybpopvalues = new double[maxnhybpopparams];
        for (int pp=0; pp<ntippopparams; pp++) {
            testtippopvalues.setParameterValue(pp, 1000+pp);
        }
        for (int pp=0; pp<nrootpopparams; pp++) {
            testrootpopvalues.setParameterValue(pp, 2000+pp);
        }
        for (int pp=0; pp<maxnhybpopparams; pp++) {
            testhybpopvalues[pp] = 3000+pp;
        }
        AlloppMulLabTree testmullabtree = new AlloppMulLabTree(adhist, tettrees, apsp,
                testtippopvalues, testrootpopvalues, testhybpopvalues);
        System.out.println(testmullabtree.asText());
        String newick = testmullabtree.mullabTreeAsNewick();
        return newick;
    }


    // for test cases
    void addSimpleNodeChildren(SimpleNode anc,  SimpleNode lch, SimpleNode rch, double minlen) {
        anc.addChild(lch);
        anc.addChild(rch);
        anc.setHeight(Math.max(lch.getHeight(), rch.getHeight()) + minlen);
    }


}
