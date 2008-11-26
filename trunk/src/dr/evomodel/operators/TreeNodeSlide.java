package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;
import jebl.util.FixedBitSet;

import java.util.Arrays;

/**
 * An operator on a species tree based on the ideas of Mau et all.
 *
 * <a href="http://citeseer.ist.psu.edu/rd/27056960%2C5592%2C1%2C0.25%2CDownload/http://citeseer.ist.psu.edu/cache/papers/cs/5768/ftp:zSzzSzftp.stat.wisc.eduzSzpubzSznewtonzSzlastfinal.pdf/mau98bayesian.pdf">
 * Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods</a>
 *
 *  @author joseph
 *         Date: 29/05/2008
 */
public class TreeNodeSlide extends SimpleMCMCOperator {
    private static final String TREE_NODE_REHEIGHT = "nodeReHeight";

    private SpeciesTreeModel tree;
    private SpeciesBindings species;

    private int[] preOrderIndexBefore;
    private int[] preOrderIndexAfter;

    private double range = 1.0;
    private boolean outgroupOnly = false;

    //private boolean verbose = false;

    public TreeNodeSlide(SpeciesTreeModel tree, SpeciesBindings species, double range, boolean outgroupOnly, double weight) {
        this.tree = tree;
        this.species = species;
        this.range = range;
        this.outgroupOnly = outgroupOnly;

        preOrderIndexBefore = new int[tree.getNodeCount()];
        Arrays.fill(preOrderIndexBefore,  -1);

        preOrderIndexAfter= new int[tree.getNodeCount()];
        Arrays.fill(preOrderIndexAfter,  -1);

        setWeight(weight);
    }

    public String getPerformanceSuggestion() {
        return "none";
    }

    public String getOperatorName() {
        return PARSER.getParserName() + "(" + tree.getId() + "," + species.getId() + ")";
    }

    public double doOperation() throws OperatorFailedException {

//            #print "operate: tree", ut.treerep(t)
     //   if( verbose)  System.out.println("  Mau at start: " + tree.getSimpleTree());

        final int count = tree.getExternalNodeCount();  assert count == species.nSpecies();

        NodeRef[] order = new NodeRef[2 * count - 1];
        boolean[] swapped = new boolean[count-1];
        mauCanonical(tree, order, swapped);


        // internal node to change
        int which = MathUtils.nextInt(count - 1);
        if( outgroupOnly ) {
           if( order[1] == tree.getRoot() ) {
               which = 0;
           } else if( order[2*count - 3] == tree.getRoot() ) {
               which = count - 2;
           }
        }
//        if( verbose)  {
//            System.out.print("order:");
//            for(int k = 0 ; k < 2*count; k += 2) {
//                System.out.print(tree.getNodeTaxon(order[k]).getId() + ((k == 2*which) ? "*" : " "));
//            }
//            System.out.println();
//        }

        FixedBitSet left = new FixedBitSet(count);
        FixedBitSet right = new FixedBitSet(count);

        for(int k = 0; k < 2*which+1; k += 2) {
           left.set(tree.speciesIndex(order[k]));
        }

        for(int k = 2*which+2; k < 2*count; k += 2) {
           right.set(tree.speciesIndex(order[k]));
        }
        
        final double limit = species.speciationUpperBound(left, right);

        double amount = (MathUtils.nextDouble() - 0.5) * limit * range;
        double h = tree.getNodeHeight(order[2*which+1]) + amount;
        if( h < 0 ) {
            h = -h;
        } else if ( h > limit ) {
            h = 2*limit - h;
        }
        final double newHeight = h; // MathUtils.nextDouble() * limit;   //h
        assert 0 < newHeight && newHeight < limit;
        
//        if( verbose)  {
//            System.out.println("limit" + limit + " newH " + newHeight);
//        }

        tree.beginTreeEdit();

        setPreorderIndices(tree, preOrderIndexBefore);

        tree.setNodeHeight(order[2*which+1], newHeight);

        mauReconstruct(tree, order, swapped);

        // restore preorder of pops -
        {
            setPreorderIndices(tree, preOrderIndexAfter);

            double[] splitPopValues = null;

            for(int k = 0; k < preOrderIndexBefore.length; ++k) {
                final int b = preOrderIndexBefore[k];
                if( b >= 0 ) {
                    final int a = preOrderIndexAfter[k];
                    if( a != b ) {
                        //if( verbose)  System.out.println("pops: " + a + " <- " + b);

                        final Parameter p1 = tree.sppSplitPopulations;
                        if( splitPopValues == null ) {
                            splitPopValues = p1.getParameterValues();
                        }
                        for(int i = 0; i < 2; ++i){
                            p1.setParameterValue(count + 2*a + i, splitPopValues[count + 2*b + i]);
                        }
                    }
                }
            }
        }

       // if( verbose)  System.out.println("  Mau after: " + tree.getSimpleTree());
       // }

        tree.endTreeEdit();

        return 0;
    }

//    static private void treeMixup(Tree tree, NodeRef node) {
//        final double h = tree.getNodeHeight(node);
//        if( ! tree.isRoot(node) )  {
//            assert tree.getBranchLength(node) ==  (tree.getNodeHeight(tree.getParent(node)) - h);
//        }
//        if( ! tree.isExternal(node) ) {
//           for(int k = 0; k < tree.getChildCount(node); ++k) {
//               assert  h >  tree.getNodeHeight(tree.getChild(node, k));
//               final NodeRef child = tree.getChild(node, k);
//               treeMixup(tree, child);
//           }
//        }
//    }

    static private void setPreorderIndices(SpeciesTreeModel tree, int[] indices) {
        setPreorderIndices(tree, tree.getRoot(), 0, indices);
    }

    static private int setPreorderIndices(SpeciesTreeModel tree, NodeRef node, int loc, int[] indices) {
        if( ! tree.isExternal(node) ) {
            int l = setPreorderIndices(tree, tree.getChild(node, 0), loc, indices);
            indices[node.getNumber()] = l;
            loc = setPreorderIndices(tree, tree.getChild(node, 1), l+1, indices);
        }
        return loc;
    }

    static private void mauCanonical(Tree tree, NodeRef[] order, boolean[] wasSwapped) {
        mauCanonicalSub(tree, tree.getRoot(), 0, order, wasSwapped);
    }

    static private int mauCanonicalSub(Tree tree, NodeRef node, int loc, NodeRef[] order, boolean[] wasSwaped) {
        if( tree.isExternal(node) ) {
            order[loc] = node;     assert (loc & 0x1) == 0;
            return loc + 1;
        }

        final boolean swap = MathUtils.nextBoolean();
        //wasSwaped[(loc-1)/2] = swap;
       
        int l = mauCanonicalSub(tree, tree.getChild(node, swap ? 1 : 0), loc, order, wasSwaped);

        order[l] = node;   assert (l & 0x1) == 1;
        wasSwaped[(l-1)/2] = swap;

        l = mauCanonicalSub(tree, tree.getChild(node, swap ? 0 : 1), l+1, order, wasSwaped);
        return l;
    }

    static private void mauReconstruct(SpeciesTreeModel tree, NodeRef[] order, boolean[] swapped) {
        final NodeRef root = mauReconstructSub(tree, 0, swapped.length, order, swapped);
        if( tree.getRoot() != root ) {
            tree.setRoot(root);
        }
    }

    static private NodeRef mauReconstructSub(SpeciesTreeModel tree, int from, int to, NodeRef[] order, boolean[] wasSwaped) {
        if( from == to ) {
            return order[2*from];
        }

        int rootIndex = -1;
        {
            double h = -1;

            for(int i = from; i < to; ++i) {
                final double v = tree.getNodeHeight(order[2 * i + 1]);
                if( h < v ) {
                    h = v;
                    rootIndex = i;
                }
            }
        }


        final NodeRef root = order[2 * rootIndex + 1];
        
        final NodeRef lchild = tree.getChild(root, 0);
        final NodeRef rchild = tree.getChild(root, 1);

        NodeRef lTargetChild = mauReconstructSub(tree, from, rootIndex, order, wasSwaped);
        NodeRef rTargetChild = mauReconstructSub(tree, rootIndex+1, to, order, wasSwaped);

        if( wasSwaped[rootIndex] ) {
            NodeRef z = lTargetChild;
            lTargetChild = rTargetChild;
            rTargetChild = z;
        }

        if( lchild != lTargetChild ) {
            tree.replaceChild(root, lchild, lTargetChild);
        }

        if( rchild != rTargetChild ) {
            tree.replaceChild(root, rchild, rTargetChild);
        }

        return root;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_NODE_REHEIGHT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            SpeciesBindings species = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
            SpeciesTreeModel tree = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);

            final double weight = xo.getDoubleAttribute("weight");
            final double range = xo.getAttribute("range", 1.0);
            if( range <= 0 || range > 1.0 ) {
                throw new XMLParseException("range out of range");
            }
            final boolean oo = xo.getAttribute("outgroup", false);
            return new TreeNodeSlide(tree, species, range, oo, weight);
        }

        public String getParserDescription() {
            return "Specialized Species tree operator, transform tree without breaking embedding of gene trees.";
        }

        public Class getReturnType() {
            return TreeNodeSlide.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule("weight"),
                    AttributeRule.newDoubleRule("range", true),
                    AttributeRule.newBooleanRule("outgroup", true),
                    
                    new ElementRule(SpeciesBindings.class),
                    new ElementRule(SpeciesTreeModel.class)
            };
        }
    };
}
