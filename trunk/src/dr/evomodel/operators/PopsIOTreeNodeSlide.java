package dr.evomodel.operators;

/*
 * User: Graham Jones
 * Date: 11/05/12
 */

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SlidableTree;
import dr.evolution.util.Taxon;
import dr.evomodel.speciation.PopsIOSpeciesBindings;
import dr.evomodel.speciation.PopsIOSpeciesTreeModel;
import dr.evomodelxml.operators.PopsIOTreeNodeSlideParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import jebl.util.FixedBitSet;


public class PopsIOTreeNodeSlide  extends SimpleMCMCOperator {
    PopsIOSpeciesTreeModel piostm;
    PopsIOSpeciesBindings piosb ;

    public PopsIOTreeNodeSlide(PopsIOSpeciesTreeModel piostm, PopsIOSpeciesBindings piosb, double weight) {
        this.piostm = piostm;
        this.piosb = piosb;
        setWeight(weight);
    }

    @Override
    public String getPerformanceSuggestion() {
        return "None";
    }

    @Override
    public String getOperatorName() {
        return PopsIOTreeNodeSlideParser.PIOTREE_NODESLIDE;
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        operateOneNodeInTree();
        return 0.0;
    }


    private int randomnode() {
        int count = piostm.getInternalNodeCount();
        int which = MathUtils.nextInt(count);
        return which;
    }





    private void operateOneNodeInTree() {
        int which = randomnode();

        NodeRef[] order = SlidableTree.Utils.mnlCanonical(piostm);

        // Find the time of the most recent gene coalescence which
        // has species to left and right of this node.
        FixedBitSet left = piosb.emptyUnion();
        FixedBitSet right = piosb.emptyUnion();
        for (int k = 0; k < 2 * which + 1; k += 2) {
            Taxon tx = piostm.getSlidableNodeTaxon(order[k]);
            left.union(piosb.tipUnionFromTaxon(tx));
        }
        for (int k = 2 * (which + 1); k < order.length; k += 2) {
            Taxon tx = piostm.getSlidableNodeTaxon(order[k]);
            right.union(piosb.tipUnionFromTaxon(tx));
        }
        double genelimit = piosb.coalescenceUpperBoundBetween(left, right);
        double newHeight = MathUtils.nextDouble() * genelimit;

        final NodeRef node = order[2 * which + 1];
        piostm.setSlidableNodeHeight(node, newHeight);
        SlidableTree.Utils.mnlReconstruct(piostm, order);
        piostm.fixupAfterNodeSlide();
    }


}
