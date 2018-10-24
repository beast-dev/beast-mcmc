package dr.evomodelxml.treelikelihood;

import com.sun.org.apache.regexp.internal.RE;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresTraitDataModel;
import dr.evomodel.treelikelihood.utilities.VarianceProportionLogger;
import dr.xml.*;

public class VarianceProportionLoggerParser extends AbstractXMLObjectParser {
    public static final String PARSER_NAME = "varianceProportionLogger";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);
        RepeatedMeasuresTraitDataModel dataModel = (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);
        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        return new VarianceProportionLogger(tree, treeLikelihood, dataModel, diffusionModel);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(RepeatedMeasuresTraitDataModel.class),
            new ElementRule(MultivariateDiffusionModel.class)
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserName() {return PARSER_NAME;}

    public String getParserDescription() {
        return "A parser to log sampling variance components";
    }

    public Class getReturnType() {
        return VarianceProportionLogger.class;
    }


}
