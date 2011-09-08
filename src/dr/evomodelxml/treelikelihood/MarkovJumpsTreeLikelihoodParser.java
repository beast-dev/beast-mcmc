package dr.evomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.xml.*;

/**
 */
public class MarkovJumpsTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "markovJumpsTreeLikelihood";

    public String getParserName() {
        return RECONSTRUCTING_TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
       throw new XMLParseException("MarkovJump functionality is only support under BEAGLE.");
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return AncestralStateTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return null;
    }
}