package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedSubstitutionModel extends BaseSubstitutionModel implements Citable {

    private List<SubstitutionModel> baseModels;
    private final int numBaseModel;
    private final int[] stateSizes;
    private final ProductChainFrequencyModel pcFreqModel;
    private final Parameter switchingRates;

    public MarkovModulatedSubstitutionModel(String name,
                                         List<SubstitutionModel> baseModels,
                                         Parameter switchingRates) {
        super(name);

        this.baseModels = baseModels;
        numBaseModel = baseModels.size();

        if (numBaseModel == 0) {
            throw new RuntimeException("May not construct ProductChainSubstitutionModel with 0 base models");
        }

        this.switchingRates = switchingRates;
        addVariable(switchingRates);

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>();
        stateSizes = new int[numBaseModel];
        stateCount = 1;
        for (int i = 0; i < numBaseModel; i++) {
            addModel(baseModels.get(i));
            freqModels.add(baseModels.get(i).getFrequencyModel());
            DataType dataType = baseModels.get(i).getDataType();
            stateSizes[i] = dataType.getStateCount();
            stateCount *= dataType.getStateCount();
        }

        pcFreqModel = new ProductChainFrequencyModel("pc",freqModels);

//        String[] codeStrings = getCharacterStrings();
//        dataType = new GeneralDataType(codeStrings);

        updateMatrix = true;

        Logger.getLogger("dr.app.beagle").info("\tConstructing a Markov-modulated Markov chain substitution model,  please cite:\n"
                + Citable.Utils.getCitationString(this));
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.OBRIEN_2009
        );
        return citations;
    }

    @Override
    protected void frequenciesChanged() {
        // Do nothing

    }

    @Override
    protected void ratesChanged() {
        updateMatrix = true;  // Lazy recompute relative rates
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == switchingRates) {
            // Update rates
            ratesChanged();
        }
        // else do nothing, action taken care of at individual base models
    }
}
