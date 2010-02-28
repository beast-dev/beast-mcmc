package dr.app.beagle.evomodel.parsers;

import dr.xml.*;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;

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

        AncestralStateBeagleTreeLikelihood[] partition = new AncestralStateBeagleTreeLikelihood[3];
        String[] labels = new String[]{FIRST, SECOND, THIRD};

        BranchRateModel testBranchRateModel = null;
        for (int i = 0; i < 3; i++) {
            partition[i] = (AncestralStateBeagleTreeLikelihood)
                    xo.getChild(labels[i]).getChild(AncestralStateBeagleTreeLikelihood.class);

            // Ensure that siteRateModel has one category
            if (partition[i].getSiteRateModel().getCategoryCount() > 1) {
                throw new XMLParseException("Robust counting currently only implemented for single category models");
            }

            // Ensure that branchRateModel is the same across all partitions
            if (testBranchRateModel == null) {
                testBranchRateModel = partition[i].getBranchRateModel();
            } else if (testBranchRateModel != partition[i].getBranchRateModel()) {
                throw new XMLParseException(
                        "Robust counting currently requires the same branch rate model for all partitions");
            }
        }

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            codons = Codons.findByName(codeStr);
        }

        String labelingString = (String) xo.getAttribute(LABELING);
        CodonLabeling codonLabeling = CodonLabeling.parseFromString(labelingString);
        if (codonLabeling == null) {
            throw new XMLParseException("Unrecognized codon labeling '" + labelingString + "'");
        }

        return new CodonPartitionedRobustCounting(
                xo.getId(),
                tree,
                partition,
                codons,
                codonLabeling);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(FIRST,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(SECOND,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(THIRD,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(TreeModel.class),
            new StringAttributeRule(GeneticCode.GENETIC_CODE,
                    "The genetic code to use",
                    GeneticCode.GENETIC_CODE_NAMES, true),
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
