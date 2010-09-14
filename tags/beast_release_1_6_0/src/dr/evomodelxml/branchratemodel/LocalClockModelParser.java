package dr.evomodelxml.branchratemodel;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.LocalClockModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class LocalClockModelParser extends AbstractXMLObjectParser {

    public static final String LOCAL_CLOCK_MODEL = "localClockModel";
    public static final String RATE = BranchRateModel.RATE;
    public static final String RELATIVE = "relative";
    public static final String CLADE = "clade";
    public static final String INCLUDE_STEM = "includeStem";
    public static final String EXCLUDE_CLADE = "excludeClade";
    public static final String EXTERNAL_BRANCHES = "externalBranches";
    public static final String BACKBONE = "backbone";

    public String getParserName() {
        return LOCAL_CLOCK_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter globalRateParameter = (Parameter) xo.getElementFirstChild(RATE);
        LocalClockModel localClockModel = new LocalClockModel(tree, globalRateParameter);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject xoc = (XMLObject) xo.getChild(i);
                if (xoc.getName().equals(CLADE)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);

                    if (taxonList.getTaxonCount() == 1) {
                        throw new XMLParseException("A local clock for a clade must be defined by at least two taxa");
                    }

                    boolean includeStem = false;
                    boolean excludeClade = false;

                    if (xoc.hasAttribute(INCLUDE_STEM)) {
                        includeStem = xoc.getBooleanAttribute(INCLUDE_STEM);
                    }

                    if (xoc.hasAttribute(EXCLUDE_CLADE)) {
                        excludeClade = xoc.getBooleanAttribute(EXCLUDE_CLADE);
                    }

                    try {
                        localClockModel.addCladeClock(taxonList, rateParameter, relative, includeStem, excludeClade);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(EXTERNAL_BRANCHES)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);


                    try {
                        localClockModel.addExternalBranchClock(taxonList, rateParameter, relative);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(BACKBONE)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);


                    try {
                        localClockModel.addBackboneClock(taxonList, rateParameter, relative);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                }

            }
        }

        System.out.println("Using local clock branch rate model.");

        return localClockModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a branch rate model that adds a delta to each terminal branch length.";
    }

    public Class getReturnType() {
        return LocalClockModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
            new ElementRule(EXTERNAL_BRANCHES,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            new ElementRule(Taxa.class, "A local clock that will be applied only to " +
                                    "the external branches for these taxa"),
                            new ElementRule(Parameter.class, "The rate parameter"),
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(CLADE,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            AttributeRule.newBooleanRule(INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel (default false)."),
                            AttributeRule.newBooleanRule(EXCLUDE_CLADE, true, "determines whether to exclude actual branches of the clade from the siteModel (default false)."),
                            new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
                            new ElementRule(Parameter.class, "The rate parameter")
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(BACKBONE,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            new ElementRule(Taxa.class, "A local clock that will be applied only to " +
                                    "the 'backbone' branches defined by these taxa"),
                            new ElementRule(Parameter.class, "The rate parameter"),
                    }, 0, Integer.MAX_VALUE),
    };
}
