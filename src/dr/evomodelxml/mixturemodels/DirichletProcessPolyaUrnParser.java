package dr.evomodelxml.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.evomodel.mixturemodels.DirichletProcessPolyaUrn;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class DirichletProcessPolyaUrnParser extends AbstractXMLObjectParser {

    public static final String DIRICHLET_PROCESS_POLYA_URN = "dirichletProcessPolyaUrn";
    public static final String BASE_MODEL = "baseModel";
    public static final String MASS = "mass";
    public static final String CATEGORIES = "categories";
    public static final String GROUP_ASSIGNMENTS = "groupAssignments";
    public static final String MAX_NUM_CAT = "maxNumberOfCategories";
    public static final String NUM_STARTING_CAT = "numStartingCategories";

    @Override
    public String getParserName() {
        return DIRICHLET_PROCESS_POLYA_URN;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
        Parameter groupAssignments = (Parameter)xo.getElementFirstChild(GROUP_ASSIGNMENTS);
        Parameter maxNumCat = (Parameter)xo.getElementFirstChild(MAX_NUM_CAT);
        int maxNumOfCat = (int)maxNumCat.getParameterValue(0);

        CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
        //ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) xo.getElementFirstChild(BASE_MODEL);
        //List<ParametricMultivariateDistributionModel> baseDistributionList = new ArrayList<>();
        //baseDistributionList.add(baseModel);
        Parameter gamma = (Parameter) xo.getElementFirstChild(MASS);
        List<Parameter> massParameterList = new ArrayList<>();
        massParameterList.add(gamma);

        List<ParametricMultivariateDistributionModel> baseDistributionList = new ArrayList<>();
        XMLObject cxo = xo.getChild(BASE_MODEL);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            ParametricMultivariateDistributionModel baseModel = (ParametricMultivariateDistributionModel) cxo.getChild(i);
            baseDistributionList.add(baseModel);
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

        CountableRealizationsParameter allParameters = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);

        return new DirichletProcessPolyaUrn(groupAssignments,categoriesParameter,
                uniquelyRealizedParameters,
                allParameters,
                baseDistributionList,
                massParameterList,
                maxNumOfCat);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{

                new ElementRule(CATEGORIES,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments

                new ElementRule(CompoundParameter.class, false), // realized parameters

                new ElementRule(MAX_NUM_CAT,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class, false)
                        }
                ),

                new ElementRule(BASE_MODEL,
                        new XMLSyntaxRule[] {
                                new ElementRule(ParametricMultivariateDistributionModel.class, 1, Integer.MAX_VALUE),
                        }
                ), // base models

                new ElementRule(MASS,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }),// gamma
                AttributeRule.newIntegerRule(NUM_STARTING_CAT, true),

        };
    }

    @Override
    public String getParserDescription() {
        return DIRICHLET_PROCESS_POLYA_URN;
    }

    @Override
    public Class getReturnType() {
        return DirichletProcessPolyaUrn.class;
    }

}