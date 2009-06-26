package dr.evomodel.continuous;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.BranchAttributeProvider;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;
import dr.math.matrixAlgebra.Vector;

/*
 * @author Marc Suchard
 */
public class BranchDirectionAttributeProvider implements BranchAttributeProvider {


    public static void main(String[] arg) {
        double[] vector = {-10,1};

        double angle = convert(vector[0],vector[1]);
      
        System.err.println("vec:   "+new Vector(vector));
        System.err.println("angle: "+angle);
    }

    private static double convert(double latValue, double longValue) {
        double angle = Math.atan2(latValue,longValue);
             if (angle < 0)
                 angle = 2*Math.PI + angle;
        return angle;
    }

    public static final String DIRECTION_PROVIDER = "branchDirections";
    public static final String FORMAT = "%5.4f";

    public BranchDirectionAttributeProvider(SampledMultivariateTraitLikelihood traitLikelihood) {
//        this.traitLikelihood = traitLikelihood;
        this.treeModel = traitLikelihood.getTreeModel();
        this.traitName = traitLikelihood.getTraitName();
        this.label = traitName + ".angle";

        double[] rootTrait = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName);
        if (rootTrait.length != 2)
            throw new RuntimeException("BranchDirectionAttributeProvider only works for 2D traits");
    }

    public String getBranchAttributeLabel() {
        return label;
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {

        if( tree != treeModel)
            throw new RuntimeException("Bad bug.");

        NodeRef parent = tree.getParent(node);
        double[] parentTrait = treeModel.getMultivariateNodeTrait(parent,traitName);
        double[] vector   = treeModel.getMultivariateNodeTrait(node,  traitName);

        vector[0] -= parentTrait[0];
        vector[1] -= parentTrait[1];

        return String.format(FORMAT,convert(vector[0],vector[1]));
    }

     public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             SampledMultivariateTraitLikelihood traitLikelihood = (SampledMultivariateTraitLikelihood)
                     xo.getChild(SampledMultivariateTraitLikelihood.class);

             return new BranchDirectionAttributeProvider(traitLikelihood);
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return new XMLSyntaxRule[] {
                 new ElementRule(SampledMultivariateTraitLikelihood.class),
             };
         }

         public String getParserDescription() {
             return null;
         }

         public Class getReturnType() {
             return NodeAttributeProvider.class;
         }

         public String getParserName() {
             return DIRECTION_PROVIDER;
         }
     };

//    private SampledMultivariateTraitLikelihood traitLikelihood;
    private TreeModel treeModel;
    private String traitName;
    private String label;


}
