package dr.app.beauti.options;

import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.generator.GeneralTraitGenerator;

import java.util.List;

/**
 *
 */
public class PartitionClockModelSubstModelLink extends PartitionOptions {
    private final BeautiOptions options;
    private final PartitionClockModel clockModel;
    private final PartitionSubstitutionModel substModel;

    public PartitionClockModelSubstModelLink(BeautiOptions options, PartitionClockModel clockModel, PartitionSubstitutionModel substModel) {
        this.options = options;
        this.clockModel = clockModel;
        this.substModel = substModel;

        initClockModelSubstModelLinkParaAndOpers();
    }

    private void initClockModelSubstModelLinkParaAndOpers() {
        // <svsGeneralSubstitutionModel idref="originModel"/>
//        createParameterAndStringOperator(OperatorType.BITFIP_IN_SUBST.toString(), getPrefix() + "trait.mu",
//                "bit Flip In Substitution Model Operator",
//                substModel.getParameter("trait.mu").getName(),  TODO trait.mu belongs Clock Model?
//                GeneralTraitGenerator.getLocationSubstModelTag(substModel), substModel.getPrefix() + substModel.getName(),
//                OperatorType.BITFIP_IN_SUBST, demoTuning, 30);

        createParameterAndStringOperator(OperatorType.BITFIP_IN_SUBST.toString(), "clock.rate",
                "bit Flip In Substitution Model Operator on clock.rate", clockModel.getParameter("clock.rate"),
                GeneralTraitGenerator.getLocationSubstModelTag(substModel), substModel.getPrefix() + substModel.getName(),
                OperatorType.BITFIP_IN_SUBST, demoTuning, 30);

    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {

    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (substModel.isActivateBSSVS()) {
            ops.add(getOperator(OperatorType.BITFIP_IN_SUBST.toString()));
        }
    }

    /////////////////////////////////////////////////////////////

    public PartitionClockModel getClockModel() {
        return clockModel;
    }

    public PartitionSubstitutionModel getSubstModel() {
        return substModel;
    }

    public String getPrefix() {
        return noDuplicatedPrefix(clockModel.getPrefix(), substModel.getPrefix());
    }

}
