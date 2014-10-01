package dr.evomodelxml.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.MSSD.SingleTipObservationProcess;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 *
 */
public class SingleTipObservationProcessParser extends AbstractXMLObjectParser {
    public final static String MODEL_NAME = "singleTipObservationProcess";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter mu = (Parameter) xo.getElementFirstChild(AnyTipObservationProcessParser.DEATH_RATE);
        Parameter lam = (Parameter) xo.getElementFirstChild(AnyTipObservationProcessParser.IMMIGRATION_RATE);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        Taxon sourceTaxon = (Taxon) xo.getChild(Taxon.class);
        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        Logger.getLogger("dr.evomodel.MSSD").info("Creating SingleTipObservationProcess model. All traits are assumed extant in " + sourceTaxon.getId() + "Initial mu = " + mu.getParameterValue(0) + " initial lam = " + lam.getParameterValue(0));

        return new SingleTipObservationProcess(treeModel, patterns, siteModel, branchRateModel, mu, lam, sourceTaxon);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the SingleTipObservationProcess for ALSTreeLikelihood calculations";
    }

    public Class getReturnType() {
        return SingleTipObservationProcess.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(PatternList.class),
            new ElementRule(Taxon.class),
            new ElementRule(SiteModel.class),
            new ElementRule(BranchRateModel.class),
            new ElementRule(AnyTipObservationProcessParser.DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(AnyTipObservationProcessParser.IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
