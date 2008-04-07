package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Package: SingleTipObservationProcess
 * Description:
 *
 *
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 2:57:14 PM
 */
public class SingleTipObservationProcess extends AnyTipObservationProcess{
    final static String MODEL_NAME="singleTipObservationProcess";
    protected Taxon sourceTaxon;

    public SingleTipObservationProcess(TreeModel treeModel, PatternList patterns, SiteModel siteModel,
                                       Parameter mu, Parameter lam, Taxon sourceTaxon) {
        super(MODEL_NAME, treeModel, patterns, siteModel, mu, lam);
        this.sourceTaxon=sourceTaxon;
    }

    public double calculateLogTreeWeight(BranchRateModel branchRateModel) {
        return -lam.getParameterValue(0)/(getAverageRate()*mu.getParameterValue(0));
    }

    /* ***************************************
     *  PARSER IMPLEMENTATION
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

    public String getParserName() { return MODEL_NAME; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter mu = (Parameter)xo.getSocketChild(DEATH_RATE);
            Parameter lam= (Parameter)xo.getSocketChild(IMMIGRATION_RATE);
            TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
            PatternList patterns = (PatternList)xo.getChild(PatternList.class);
            Taxon sourceTaxon = (Taxon)xo.getChild(Taxon.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
            Logger.getLogger("dr.treelikelihood").info("Creating SingleTipObservationProcess model. All traits are assumed extant in "+sourceTaxon.getId()+"Initial mu = " + mu.getParameterValue(0)+" initial lam = " + lam.getParameterValue(0));

            return new SingleTipObservationProcess(treeModel, patterns, siteModel, mu, lam, sourceTaxon);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an instance of the SingleTipObservationProcess for ALSTreeLikelihood calculations";
        }

        public Class getReturnType() { return SingleTipObservationProcess.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(TreeModel.class),
            new ElementRule(PatternList.class),
            new ElementRule(Taxon.class),
            new ElementRule(SiteModel.class),
            new ElementRule(DEATH_RATE, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
            new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
        };

    };
}
