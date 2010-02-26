package dr.app.beagle.evomodel.parsers;

import dr.xml.*;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard 
 */
public class CodonPartitionedRobustCountingParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "codonPartitionedRobustCounting";
    public static final String FIRST = "firstPosition";
    public static final String SECOND = "secondPosition";
    public static final String THIRD = "thirdPosition";
    public static final String LABELING = "labeling";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
   
        AncestralStateBeagleTreeLikelihood partition0 = (AncestralStateBeagleTreeLikelihood)
                xo.getChild(FIRST).getChild(AncestralStateBeagleTreeLikelihood.class);
        AncestralStateBeagleTreeLikelihood partition1 = (AncestralStateBeagleTreeLikelihood)
                xo.getChild(SECOND).getChild(AncestralStateBeagleTreeLikelihood.class);
        AncestralStateBeagleTreeLikelihood partition2 = (AncestralStateBeagleTreeLikelihood)
                xo.getChild(THIRD).getChild(AncestralStateBeagleTreeLikelihood.class);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        String labelingString = (String) xo.getAttribute(LABELING);
        CodonLabeling codonLabeling = CodonLabeling.parseFromString(labelingString);
        if (codonLabeling == null) {
            throw new XMLParseException("Unrecognized codon labeling '" + labelingString +"'");
        }

        return new CodonPartitionedRobustCounting(xo.getId(), tree, partition0,  partition1,  partition2, codonLabeling);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(FIRST,
                    new XMLSyntaxRule[] {
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(SECOND,
                    new XMLSyntaxRule[] {
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(THIRD,
                    new XMLSyntaxRule[] {
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(TreeModel.class),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "A parser to specify robust counting procedures on codon partitioned models";
    }

    public Class getReturnType() {
        return CodonPartitionedRobustCounting.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
