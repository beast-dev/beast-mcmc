package dr.app.beauti.components.continuous;

import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.evolution.datatype.ContinuousDataType;
import no.uib.cipr.matrix.SymmDenseEVD;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */

public class ContinuousComponentOptions implements ComponentOptions {

    public final static String PRECISION_GIBBS_OPERATOR = "precisionGibbsOperator";
    public final static String HALF_DF = "halfDF";
    public final static String STDEV = "stdev";
    public final static String LAMBDA = "lambda";
    public final static String RRW_CATEGORIES = "rrwCategories";

    final private BeautiOptions options;

	public ContinuousComponentOptions(final BeautiOptions options) {
		this.options = options;
	}

	public void createParameters(ModelOptions modelOptions) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            String prefix = partitionData.getName() + ".";

            if (!modelOptions.parameterExists(prefix + HALF_DF)) {
                modelOptions.createParameterGammaPrior(prefix + HALF_DF, "half DF of 1 parameter gamma distributed RRW",
                        PriorScaleType.NONE, 0.5, 0.001, 1000.0, false);
                modelOptions.createScaleOperator(prefix + HALF_DF, modelOptions.demoTuning, 1.0);
            }

            if (!modelOptions.parameterExists(prefix + STDEV)) {
                modelOptions.createParameterExponentialPrior(prefix + STDEV, "standard deviation of lognormal distributed RRW",
                        PriorScaleType.NONE, 1.0 / 3.0, 1.0 / 3.0, 0.0);
                modelOptions.createScaleOperator(prefix + STDEV, modelOptions.demoTuning, 5.0);
            }

            if (!modelOptions.parameterExists(prefix + LAMBDA)) {
                modelOptions.createParameterBetaDistributionPrior(prefix + LAMBDA,
                        "phylogenetic signal parameter",
                        0.5, 2.0, 2.0, 0.0);
                modelOptions.createOperator(prefix + LAMBDA, OperatorType.RANDOM_WALK_ABSORBING, 0.3, 10.0);
            }

            if (!modelOptions.parameterExists(prefix + "swap." + RRW_CATEGORIES)) {
                modelOptions.createParameter(prefix + RRW_CATEGORIES, "relaxed random walk rate categories");

                modelOptions.createOperator(prefix + "swap." + RRW_CATEGORIES, prefix + RRW_CATEGORIES, "Performs a swap of RRW rate categories",
                        prefix + RRW_CATEGORIES, OperatorType.SWAP, 1, 30);
                modelOptions.createOperator(prefix + "draw." + RRW_CATEGORIES, prefix + RRW_CATEGORIES, "Performs an integer uniform draw of RRW rate categories",
                        prefix + RRW_CATEGORIES, OperatorType.INTEGER_UNIFORM, 1, 10);
            }

        }
	}

	public void selectOperators(ModelOptions modelOptions, List<Operator> ops) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            String prefix = partitionData.getName() + ".";

            if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.GAMMA_RRW) {
                ops.add(modelOptions.getOperator(prefix + HALF_DF));
            } else if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
                ops.add(modelOptions.getOperator(prefix + STDEV));
            }
            if (useLambda(partitionData.getPartitionSubstitutionModel())) {
                ops.add(modelOptions.getOperator(prefix + LAMBDA));
            }
            if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() != ContinuousSubstModelType.HOMOGENOUS) {
                ops.add(modelOptions.getOperator(prefix + "swap." + RRW_CATEGORIES));
                ops.add(modelOptions.getOperator(prefix + "draw." + RRW_CATEGORIES));
            }
        }
	}

	public void selectParameters(ModelOptions modelOptions, List<Parameter> params) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(ContinuousDataType.INSTANCE)) {
            if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.GAMMA_RRW) {
                params.add(modelOptions.getParameter(partitionData.getName() + "." + HALF_DF));
            } else if (partitionData.getPartitionSubstitutionModel().getContinuousSubstModelType() == ContinuousSubstModelType.LOGNORMAL_RRW) {
                params.add(modelOptions.getParameter(partitionData.getName() + "." + STDEV));
            }
            if (useLambda(partitionData.getPartitionSubstitutionModel())) {
                params.add(modelOptions.getParameter(partitionData.getName() + "." + LAMBDA));
            }
        }
	}

	public void selectStatistics(ModelOptions modelOptions,
			List<Parameter> stats) {
		// Do nothing
	}

	public BeautiOptions getOptions() {
		return options;
	}

    public boolean useLambda(PartitionSubstitutionModel model) {
        Boolean useLambda = useLambdaMap.get(model);
        if (useLambda != null) {
            return useLambda;
        }
        return false;
    }

    public void setUseLambda(PartitionSubstitutionModel model, boolean useLambda) {
        useLambdaMap.put(model, useLambda);
    }

    final private Map<PartitionSubstitutionModel, Boolean> useLambdaMap = new HashMap<PartitionSubstitutionModel, Boolean>();
}