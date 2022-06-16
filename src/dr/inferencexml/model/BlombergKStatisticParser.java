package dr.inferencexml.model;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.BlombergKStatistic;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TRAIT_NAME;

public class BlombergKStatisticParser extends AbstractXMLObjectParser {
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        String traitName = xo.getStringAttribute(TRAIT_NAME);
        return new BlombergKStatistic(treeDataLikelihood, traitName);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TRAIT_NAME)
        };
    }

    @Override
    public String getParserDescription() {
        return "Blomberg's K statistic of phylogenetic signal";
    }

    @Override
    public Class getReturnType() {
        return BlombergKStatistic.class;
    }

    @Override
    public String getParserName() {
        return BlombergKStatistic.BLOMBERGS_K;
    }
}
