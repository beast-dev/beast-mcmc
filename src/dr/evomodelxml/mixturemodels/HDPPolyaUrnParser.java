package dr.evomodelxml.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.evomodel.mixturemodels.HDPPolyaUrn;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class HDPPolyaUrnParser extends AbstractXMLObjectParser {

    public static final String HDP_POLYA_URN = "hdpPolyaUrn";
    public static final String COMMON_BASE_DIST_BASE_DIST = "commonBaseDistBaseDist";
    public static final String COMMON_BASE_DIST_MASS = "commonBaseDistMass";
    public static final String COMMON_MASS = "commonMass";
    public static final String CATEGORIES = "categories";
    public static final String GROUP_ASSIGNMENTS = "groupAssignments";
    public static final String MAX_NUM_CAT = "maxNumberOfCategories";
    public static final String MAX_GROUPS = "maxGroups";
    public static final String NUM_STARTING_CAT = "numStartingCategories";
    public static final String URV = "uniquelyRealizedValues";
    public static final String TABLE_COUNTS = "tableCounts";
    //public static final String COMMON_BASE_DIST_WEIGHTS = "commonBaseDistWeights";
    public static final String STICK_PROPORTIONS = "stickProportions";
    public static final String IS_IHMM = "isIHMM";

    @Override
    public String getParserName() {
        return HDP_POLYA_URN;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
        Parameter groupAssignments = (Parameter)xo.getElementFirstChild(GROUP_ASSIGNMENTS);
        Parameter maxNumCat = (Parameter)xo.getElementFirstChild(MAX_NUM_CAT);
        int maxNumOfCat = (int)maxNumCat.getParameterValue(0);
        int maxGroupCount = maxNumOfCat + 1;
        CompoundParameter tableCounts = (CompoundParameter)xo.getElementFirstChild(TABLE_COUNTS);
        //CompoundParameter commonBaseDistWeights = (CompoundParameter)xo.getElementFirstChild(COMMON_BASE_DIST_WEIGHTS);
        CompoundParameter stickProportions = (CompoundParameter)xo.getElementFirstChild(STICK_PROPORTIONS);
        CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getElementFirstChild(URV);
        Parameter commonBaseDistMass = (Parameter)xo.getElementFirstChild(COMMON_BASE_DIST_MASS);
        //ParametricMultivariateDistributionModel commonBaseDistBaseDist = (ParametricMultivariateDistributionModel) xo.getElementFirstChild(COMMON_BASE_DIST_BASE_DIST);
        Parameter commonMass = (Parameter) xo.getElementFirstChild(COMMON_MASS);
        List<Parameter> massParameterList = new ArrayList<>();
        massParameterList.add(commonMass);

        List<ParametricMultivariateDistributionModel> commonBaseDistBaseDist = new ArrayList<>();
        XMLObject cxo = xo.getChild(COMMON_BASE_DIST_BASE_DIST);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) cxo.getChild(i);
            commonBaseDistBaseDist.add(baseModel);
        }

        if(stickProportions.getParameterCount() != maxNumOfCat){
            throw new XMLParseException("stickProportions must have as many parameters as max number of categories");
        }

        int numStartingCategories = 0;
        if(xo.hasAttribute(NUM_STARTING_CAT)) {
            numStartingCategories = xo.getIntegerAttribute(NUM_STARTING_CAT);
            double numCat = numStartingCategories;
            if(numStartingCategories > 0){
                double[] startingCatProbs = new double[numStartingCategories];
                for(int k = 0; k < numStartingCategories; k++){
                    startingCatProbs[k] = 1/numCat;
                }
                for(int i = 0; i < categoriesParameter.getSize(); i++) {
                    categoriesParameter.setParameterValue(i, MathUtils.randomChoicePDF(startingCatProbs));
                }
            }

        }

        boolean isIHMM = true;
        if(xo.hasAttribute(IS_IHMM)) {
            isIHMM = xo.getBooleanAttribute(IS_IHMM);
        }

        if(xo.hasAttribute(MAX_GROUPS)){
            maxGroupCount = xo.getIntegerAttribute(MAX_GROUPS);
        }

        if(tableCounts.getParameterCount() != maxGroupCount){
            throw new XMLParseException("tableCounts must have number of parameters equal to max number of groups");
        }

        CountableRealizationsParameter allParameters = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);

        return new HDPPolyaUrn(groupAssignments,categoriesParameter,
                uniquelyRealizedParameters,
                allParameters,
                massParameterList,
                commonBaseDistBaseDist,
                commonBaseDistMass,
                maxGroupCount,
                maxNumOfCat,
                isIHMM,
                tableCounts,
                stickProportions);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{

                new ElementRule(CATEGORIES,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments
                new ElementRule(MAX_NUM_CAT,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class, false)
                        }
                ),

                //new ElementRule(BASE_MODEL,
                //        new XMLSyntaxRule[] {
                //                new ElementRule(ParametricMultivariateDistributionModel.class, 1, Integer.MAX_VALUE),
                //        }
                //), // base models
                new ElementRule(GROUP_ASSIGNMENTS,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
                new ElementRule(URV,
                        new XMLSyntaxRule[] { new ElementRule(CompoundParameter.class, false) }),
                new ElementRule(TABLE_COUNTS,
                        new XMLSyntaxRule[] { new ElementRule(CompoundParameter.class, false) }),
                new ElementRule(STICK_PROPORTIONS,
                        new XMLSyntaxRule[] { new ElementRule(CompoundParameter.class, false) }),
                new ElementRule(COMMON_BASE_DIST_MASS,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
                new ElementRule(COMMON_MASS,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),
                AttributeRule.newIntegerRule(NUM_STARTING_CAT, true),
                AttributeRule.newBooleanRule(IS_IHMM, true),
                AttributeRule.newIntegerRule(MAX_GROUPS, true)
        };
    }

    @Override
    public String getParserDescription() {
        return HDP_POLYA_URN;
    }

    @Override
    public Class getReturnType() {
        return HDPPolyaUrn.class;
    }

}