package dr.evomodelxml.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.evomodel.mixturemodels.HMMPolyaUrn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class HMMPolyaUrnParser extends AbstractXMLObjectParser {

    public static final String HMM_POLYA_URN = "hmmPolyaUrn";
    //public static final String BASE_MODEL = "baseModel";
    //public static final String MASS = "mass";
    public static final String CATEGORIES = "categories";
    public static final String GROUP_ASSIGNMENTS = "groupAssignments";
    public static final String NUM_CAT = "numberOfCategories";
    public static final String DIRICHLET_PRIOR_CONCENTRATIONS = "dirichletPriorConcentrations";

    @Override
    public String getParserName() {
        return HMM_POLYA_URN;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter categoriesParameter =  (Parameter)xo.getElementFirstChild(CATEGORIES);
        Parameter groupAssignments = (Parameter)xo.getElementFirstChild(GROUP_ASSIGNMENTS);
        Parameter numCat = (Parameter)xo.getElementFirstChild(NUM_CAT);
        int numOfCat = (int)numCat.getParameterValue(0);
        int numGroups = 1 + numOfCat;

        CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);

        XMLObject cxo = xo.getChild(DIRICHLET_PRIOR_CONCENTRATIONS);

        List<Parameter> dpConcentrations = new ArrayList<Parameter>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Object testObject = cxo.getChild(i);
            if (testObject instanceof Parameter) {
                Parameter dpc = (Parameter) testObject;
                dpConcentrations.add(dpc);
                if(dpc.getDimension() != numOfCat){
                    throw new XMLParseException("Each dirichletPriorConcentration parameter must have dimension equal to number of categories");
                }
            }
        }

        if(dpConcentrations.size() != numGroups){
            throw new XMLParseException("Must have same number of groups and dirichletPriorConcentration parameters");
        }

        CountableRealizationsParameter allParameters = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);

        return new HMMPolyaUrn(groupAssignments,categoriesParameter,
                uniquelyRealizedParameters,
                allParameters,
                dpConcentrations,
                numGroups,
                numOfCat
        );
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{

                new ElementRule(CATEGORIES,
                        new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }), // categories assignments
                new ElementRule(CompoundParameter.class, false), // realized parameters
                new ElementRule(NUM_CAT,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class, false)}),
                new ElementRule(DIRICHLET_PRIOR_CONCENTRATIONS, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)
                }),

        };
    }

    @Override
    public String getParserDescription() {
        return HMM_POLYA_URN;
    }

    @Override
    public Class getReturnType() {
        return HMMPolyaUrn.class;
    }

}
