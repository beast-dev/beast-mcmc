package dr.evomodel.mixturemodels;

import dr.evomodel.branchmodel.lineagespecific.CountableRealizationsParameter;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.distribution.PointMassMixtureDistributionModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import java.util.ArrayList;
import java.util.List;

public class HMMPolyaUrn extends GenPolyaUrnProcessPrior {

    public HMMPolyaUrn(Parameter groupAssignments,
                       Parameter categoriesParameter,
                       CompoundParameter uniquelyRealizedParameters,
                       CountableRealizationsParameter allParameters,
                       List<Parameter> dirichletPriorConcentrations,
                       int numberOfGroups,
                       int numberOfCats
    ) {
        super(groupAssignments,categoriesParameter,uniquelyRealizedParameters, allParameters, null,
                null,dirichletPriorConcentrations, null,
                null,numberOfGroups,numberOfCats, null, null,false, true);
    }


    // Appropriate for Dirichlet Process, override for other processes
    public List<ParametricMultivariateDistributionModel> setUpBaseDistributionList(){
        List<ParametricMultivariateDistributionModel> baseDistList = new ArrayList<>();
        // K categories with indices 0,1,...,K-1
        // K+1 groups with indices 0,1,...,K

        // first handle F_0, which does not have weights coming from Dirichlet prior
        //Parameter initialWeights = new Parameter.Default(maxCategoryCount);
        //for(int k = 0; k < maxCategoryCount; k++){
        //    initialWeights.setParameterValue(k,1/maxCategoryCount);
        //}

        //PointMassMixtureDistributionModel initialBaseDist = new PointMassMixtureDistributionModel(initialWeights, uniquelyRealizedParameters, true);
        //baseDistList.add(initialBaseDist);

        for(int i = 0; i < maxGroupCount; i++){
            PointMassMixtureDistributionModel baseDist = new PointMassMixtureDistributionModel(dpConcentrations.get(i),uniquelyRealizedParameters,false);
            baseDistList.add(baseDist);
        }
        return baseDistList;
    }

    // Appropriate for Dirichlet Process, override for other processes
    public List<Parameter> setUpMassParameterList(){
        List<Parameter> mpList = new ArrayList<>();
        // M_0 can be arbitrary. We set it to 1
        Parameter mp = new Parameter.Default(1.0);
        mpList.add(mp);

        for(int i = 1; i < maxGroupCount; i++){
            double sum = 0;
            for(int j = 0; j < maxCategoryCount; j++){
                sum = sum + dpConcentrations.get(i).getParameterValue(j);
            }
            Parameter massParam = new Parameter.Default(sum);
            mpList.add(massParam);
        }
        return mpList;
    }

    // N sites with indices 0,1,...,N-1
    // K categories with indices 0,1,...,K-1
    // K+1 groups with indices 0,1,..,K
    // g_0 == 0
    // For i > 0, g_i is 1 + the category associated with site i-1
    // if site i-1 has category k, g_i == k+1
    /*
    public Parameter computeGroupAssignments(Parameter catParam){
        Parameter groupAssign = new Parameter.Default(catParam.getSize());
        groupAssign.setParameterValue(0,0.0);
        for(int i = 1; i < groupAssign.getSize(); i++){
            groupAssign.setParameterValue(i,1.0 + catParam.getParameterValue(i-1));
        }
        return groupAssign;
    }
    */

    public void computeGroupAssignments(Parameter catParam){
        //Parameter groupAssign = new Parameter.Default(catParam.getSize());
        //groupAssign.setParameterValue(0,0.0);
        for(int i = 1; i < groupAssignments.getSize(); i++){
            groupAssignments.setParameterValue(i,1.0 + catParam.getParameterValue(i-1));
        }
    }

    public void setGroup(Parameter gAssignments, int indexToSet, int catPrecedingSite){
        gAssignments.setParameterValue(indexToSet, 1 + catPrecedingSite);
    }

    public int getGroupAssignment(int catPrecedingSite){
        return 1 + catPrecedingSite;
    }

    public int getOrder(){
        return 1;
    }

    public boolean isFinite() {
        return true;
    }

    public boolean isHMM() {
        return true;
    }


    public double getSingletonProbability(Parameter param, List<ParametricMultivariateDistributionModel> baseDistList, int groupNum){
        //return getLogDensity(param, groupNum);
        double val[] = param.getParameterValues();
        return Math.exp(baseDistList.get(groupNum).logPdf(val));
    }

    public double getLogGammaRatio(int eta, int groupNum, int catNum){
        double returnVal = 0;
        for(int k = 1; k < eta; k++){
            Parameter p = uniquelyRealizedParameters.getParameter(catNum);
            returnVal = returnVal + Math.log(getMassParam(groupNum)*getSingletonProbability(p,baseDistributionList,groupNum) + k);
        }
        return returnVal;
    }

    public double getLogDensity(Parameter parameter, int groupNumber) {
        double value[] = parameter.getParameterValues();
        return baseDistributionList.get(groupNumber).logPdf(value);
    }

    @Override
    public Model getModel() {
        return this;
    }

}
