package dr.evomodel.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import java.util.List;

public class DirichletProcessPolyaUrn extends GenPolyaUrnProcessPrior {

    public DirichletProcessPolyaUrn(Parameter groupAssignments,
                                    Parameter categoriesParameter,
                                    CompoundParameter uniquelyRealizedParameters,
                                    CountableRealizationsParameter allParameters,
                                    List<ParametricMultivariateDistributionModel> baseDistributionList,
                                    List<Parameter> massParameterList,
                                    int maxCat
    ) {
        super(groupAssignments,categoriesParameter,uniquelyRealizedParameters,allParameters, baseDistributionList,massParameterList, null, null, null, 1,maxCat,null, null,false, false);
    }


}
