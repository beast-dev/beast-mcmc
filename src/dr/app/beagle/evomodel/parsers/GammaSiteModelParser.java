package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class GammaSiteModelParser extends AbstractXMLObjectParser {

    public static final String SITE_MODEL = SiteModel.SITE_MODEL;
    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String BRANCH_SUBSTITUTION_MODEL = "branchSubstitutionModel";
    public static final String SUBSTITUTION_RATE = "mutationRate";
    public static final String RELATIVE_RATE = "relativeRate";
    public static final String GAMMA_SHAPE = "gammaShape";
    public static final String GAMMA_CATEGORIES = "gammaCategories";
    public static final String PROPORTION_INVARIANT = "proportionInvariant";

    public String getParserName() {
        return SITE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		String msg = "";
		SubstitutionModel substitutionModel = null;

		boolean CHECK_BRANCH_SUBSTITUTION_MODEL = false;

		for (int i = 0; i < xo.getChildCount(); i++) {

			XMLObject cxo = (XMLObject) xo.getChild(i);

			if (cxo instanceof BranchSubstitutionModel) {

				CHECK_BRANCH_SUBSTITUTION_MODEL = true;

			}// END: BSM check
		}// END: children loop

        Parameter muParam = null;
        if (xo.hasChildNamed(SUBSTITUTION_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATE);

            msg += "\n  with initial substitution rate = " + muParam.getParameterValue(0);
        } else if (xo.hasChildNamed(RELATIVE_RATE)) {
            muParam = (Parameter) xo.getElementFirstChild(RELATIVE_RATE);

            msg += "\n  with initial relative rate = " + muParam.getParameterValue(0);
        }

        Parameter shapeParam = null;
        int catCount = 4;
        if (xo.hasChildNamed(GAMMA_SHAPE)) {
            XMLObject cxo = xo.getChild(GAMMA_SHAPE);
            catCount = cxo.getIntegerAttribute(GAMMA_CATEGORIES);
            shapeParam = (Parameter) cxo.getChild(Parameter.class);

            msg += "\n  " + catCount + " category discrete gamma with initial shape = " + shapeParam.getParameterValue(0);
        }

        Parameter invarParam = null;
        if (xo.hasChildNamed(PROPORTION_INVARIANT)) {
            invarParam = (Parameter) xo.getElementFirstChild(PROPORTION_INVARIANT);
            msg += "\n  initial proportion of invariant sites = " + invarParam.getParameterValue(0);
        }

        if (msg.length() > 0) {
            Logger.getLogger("dr.evomodel").info("Creating site model: " + msg);
        } else {
            Logger.getLogger("dr.evomodel").info("Creating site model.");
        }

        GammaSiteRateModel siteRateModel = new GammaSiteRateModel(SITE_MODEL, muParam, shapeParam, catCount, invarParam);

        if(!CHECK_BRANCH_SUBSTITUTION_MODEL) {
        // set this to pass it along to the TreeLikelihoodParser...
        siteRateModel.setSubstitutionModel(substitutionModel);
        
        }
        
        return siteRateModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SiteModel that has a gamma distributed rates across sites";
    }

    public Class<GammaSiteRateModel> getReturnType() {
        return GammaSiteRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
    		
          new XORRule(
          new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                  new ElementRule(SubstitutionModel.class)
                  }),
          new ElementRule(BRANCH_SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                  new ElementRule(BranchSubstitutionModel.class)
                  }), false
           ),
    		
            new XORRule(
                    new ElementRule(SUBSTITUTION_RATE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(RELATIVE_RATE, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }), true
            ),
            
            new ElementRule(GAMMA_SHAPE, new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule(GAMMA_CATEGORIES, true),
                    new ElementRule(Parameter.class)
            }, true),
            
            new ElementRule(PROPORTION_INVARIANT, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true)
            
    };
    
}//END: class
