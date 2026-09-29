package dr.evomodelxml.mixturemodels;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.mixturemodels.DataSquashingOperator;
import dr.evomodel.mixturemodels.GenPolyaUrnProcessPrior;
import dr.inference.model.CompoundLikelihood;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class DataSquashingOperatorParser extends AbstractXMLObjectParser {

    public static final String DATA_SQUASHING_OPERATOR= "dataSquashingOperator";
    public static final String DATA_LOG_LIKELIHOOD = "dataLogLikelihood";
    public static final String GROUP_ASSIGNMENTS = "groupAssignments";
    public static final String MH_STEPS = "mhSteps";
    public static final String CATEGORIES = "categories";
    public static final String CYCLICAL = "cyclical";
    public static final String DIST_METHOD = "distMethod";
    public static final String FIXED_NUMBER = "fixedNumber";
    public static final String STRICT_CUTOFF = "strictCutoff";
    public static final String EPSILON = "epsilon";
    public static final String SAMPLE_PROPORTION = "sampleProportion";
    public static final String MAX_NEW_CAT = "maxNewCat";
    public static final String SITE_MODELS = "siteModels";
    public static final String OLD = "old";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean old = xo.getAttribute(OLD, false);

        GenPolyaUrnProcessPrior gpuprocess = (GenPolyaUrnProcessPrior) xo.getChild(GenPolyaUrnProcessPrior.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);

        CompoundLikelihood cl = null;
        TreeDataLikelihood tdl = null;

        if(old){
            tdl = (TreeDataLikelihood) xo.getElementFirstChild(DATA_LOG_LIKELIHOOD);
        }else{
            cl = (CompoundLikelihood) xo.getElementFirstChild(DATA_LOG_LIKELIHOOD);
        }

        boolean cyclical = xo.getBooleanAttribute(CYCLICAL);

        int distMethod = xo.getIntegerAttribute(DIST_METHOD);

        if(distMethod != 1 && distMethod != 2){
            throw new XMLParseException("Must specify valid distMethod: 1 or 2");
        }

        int fixedNumber = 0;
        if(xo.hasAttribute(FIXED_NUMBER)){
            fixedNumber = xo.getIntegerAttribute(FIXED_NUMBER);
        }

        boolean strictCutoff = false;
        if(xo.hasAttribute(STRICT_CUTOFF)) {
            strictCutoff = xo.getBooleanAttribute(STRICT_CUTOFF);
        }

        //if(fixedNumber > 0){
        //    strictCutoff = true;
        //}

        double epsilon = xo.getDoubleAttribute(EPSILON);

        if(epsilon >= 1 || epsilon <= 0){
            throw new XMLParseException("epsilon must be between grater than 0 and less than or equal to 1");
        }

        double sampleProportion = 0.25;
        if(xo.hasAttribute(SAMPLE_PROPORTION)) {
            sampleProportion = xo.getDoubleAttribute(SAMPLE_PROPORTION);
        }

        int maxNewCat = 1;
        if(xo.hasAttribute(MAX_NEW_CAT)){
            maxNewCat = xo.getIntegerAttribute(MAX_NEW_CAT);
        }

        if(sampleProportion >= 1 || sampleProportion <= 0){
            throw new XMLParseException("sampleProportion must be between grater than 0 and less than or equal to 1");
        }

        int M = xo.getIntegerAttribute(MH_STEPS);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Tree treeModel = (Tree) xo.getChild(Tree.class);

        List<SiteRateModel> siteRateModelList = new ArrayList<SiteRateModel>();

        XMLObject cxo = xo.getChild(SITE_MODELS);

        if(cxo != null) {
            for (int j = 0; j < cxo.getChildCount(); j++) {
                Object testObject = cxo.getChild(j);
                if (testObject instanceof GammaSiteRateModel) {
                    GammaSiteRateModel srm = (GammaSiteRateModel) testObject;
                    siteRateModelList.add(srm);
                }
            }
        }

        return new DataSquashingOperator(gpuprocess,
                tdl,
                cl,
                siteRateModelList,
                treeModel,
                patternList,
                M,
                weight,
                cyclical,
                distMethod,
                epsilon,
                sampleProportion,
                fixedNumber,
                strictCutoff,
                maxNewCat,
                old
        );

    }// END: parseXMLObject

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(GenPolyaUrnProcessPrior.class, false),
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newIntegerRule(DIST_METHOD, false),
                AttributeRule.newBooleanRule(STRICT_CUTOFF, true),
                AttributeRule.newBooleanRule(CYCLICAL, false),
                AttributeRule.newDoubleRule(EPSILON, false),
                AttributeRule.newDoubleRule(SAMPLE_PROPORTION, true),
                AttributeRule.newIntegerRule(FIXED_NUMBER, true),
                AttributeRule.newIntegerRule(MAX_NEW_CAT, true)
        };

    }// END: getSyntaxRules

    @Override
    public String getParserName() {
        return DATA_SQUASHING_OPERATOR;
    }

    @Override
    public String getParserDescription() {
        return DATA_SQUASHING_OPERATOR;
    }

    @Override
    public Class getReturnType() {
        return DataSquashingOperator.class;
    }

}// END: class