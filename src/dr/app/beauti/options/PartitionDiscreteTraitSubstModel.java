package dr.app.beauti.options;

import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.operators.RateBitExchangeOperator;

import java.util.List;

/**
 * @author Walter Xie
 */
public class PartitionDiscreteTraitSubstModel extends PartitionSubstitutionModel {
    public static enum LocationSubstModelType {
        SYM_SUBST("Symmetric substitution model"),
        ASYM_SUBST("Asymmetric substitution model");

        LocationSubstModelType(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

//    public static final String PHYLOGEOGRAPHIC = "phylogeographic ";

    private LocationSubstModelType locationSubstType = LocationSubstModelType.SYM_SUBST;
    private boolean activateBSSVS = false;

    public PartitionDiscreteTraitSubstModel(BeautiOptions options, TraitData partition) {
        super(options, partition.getName(), new GeneralDataType());
    }

    public PartitionDiscreteTraitSubstModel(BeautiOptions options, String name, PartitionDiscreteTraitSubstModel source) {
        super(options, name, source);
    }

    public PartitionDiscreteTraitSubstModel(BeautiOptions options, String name, DataType dataType) {
        super(options, name, dataType);
    }


    protected void initSubstModelParaAndOpers() {

        createParameterUniformPrior("frequencies", getName() + " base frequencies", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createCachedGammaPrior("rates", "location substitution model rates",
                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0, 1.0, 0, Double.POSITIVE_INFINITY, false);
        createParameter("indicators", "location substitution model rate indicators");
//
//        createParameterExponentialPrior("mu", PHYLOGEOGRAPHIC + "mutation rate parameter",
//                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.1, 1.0, 0.0, 0.0, 10.0);

        createDiscreteStatistic("nonZeroRates", "for mutation rate parameter (if BSSVS was selected)");  // BSSVS was selected

        createOperator("rates", OperatorType.SCALE_INDEPENDENTLY, demoTuning, 30);
//        createScaleOperator("mu", "for mutation rate parameter", demoTuning, 10);

        createOperator("indicators", OperatorType.BITFLIP, -1.0, 30);
// TODO       createTagInsideOperator(OperatorType.BITFIP_IN_SUBST.toString(), "mu",
//                "bit Flip In Substitution Model Operator", "mu", OperatorType.BITFIP_IN_SUBST,
//                SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL, PREFIX_ + AbstractSubstitutionModel.MODEL, demoTuning, 30);
        // <svsGeneralSubstitutionModel idref="originModel"/>
        createOperatorUsing2Parameters(RateBitExchangeOperator.OPERATOR_NAME, "(indicators, rates)",
                "rateBitExchangeOperator (If both BSSVS and asymmetric subst selected)",
                "indicators", "rates", OperatorType.RATE_BIT_EXCHANGE, -1.0, 6.0);
    }


    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    @Override
    public void selectParameters(List<Parameter> params) {
        params.add(getParameter("frequencies"));
        params.add(getParameter("rates"));
//        params.add(getParameter("mu"));

        if (activateBSSVS) {
            params.add(getParameter("indicators"));
            params.add(getParameter("nonZeroRates"));
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    @Override
    public void selectOperators(List<Operator> ops) {
        ops.add(getOperator("rates"));
//        ops.add(getOperator("mu"));

        if (activateBSSVS) {
            ops.add(getOperator("indicators"));
            ops.add(getOperator(OperatorType.BITFIP_IN_SUBST.toString()));

            if (locationSubstType == LocationSubstModelType.ASYM_SUBST)
                ops.add(getOperator(RateBitExchangeOperator.OPERATOR_NAME));
        }

    }


    /////////////////////////////////////////////////////////////

    public LocationSubstModelType getLocationSubstType() {
        return locationSubstType;
    }

    public void setLocationSubstType(LocationSubstModelType locationSubstType) {
        this.locationSubstType = locationSubstType;
    }

    public boolean isActivateBSSVS() {
        return activateBSSVS;
    }

    public void setActivateBSSVS(boolean activateBSSVS) {
        this.activateBSSVS = activateBSSVS;
    }


}
