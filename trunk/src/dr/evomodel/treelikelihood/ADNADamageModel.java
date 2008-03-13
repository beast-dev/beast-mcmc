package dr.evomodel.treelikelihood;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.PatternList;
import dr.evolution.util.TaxonList;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class ADNADamageModel extends TipPartialsModel {

    public static final String ADNA_DAMAGE_MODEL = "aDNADamageModel";
    public static final String BASE_DAMAGE_RATE = "baseDamageRate";
    public static final String AGE_RATE_FACTOR = "ageRateFactor";

    public ADNADamageModel(TreeModel treeModel, PatternList patternList, Parameter baseDamageRateParameter, Parameter ageRateFactorParameter) {
        super(ADNA_DAMAGE_MODEL);

        this.treeModel = treeModel;
        this.baseDamageRateParameter = baseDamageRateParameter;
        this.ageRateFactorParameter = ageRateFactorParameter;

        DataType dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

        int extNodeCount = treeModel.getExternalNodeCount();

        states = new int[extNodeCount][];
        partials = new double[extNodeCount][];

        try {

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                }

                setStates(patternList, index, i);
            }
        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
    }

    /**
     * Sets the partials from a sequence in an alignment.
     * @param patternList
     * @param sequenceIndex
     * @param nodeIndex
     */
    private final void setStates(PatternList patternList, int sequenceIndex, int nodeIndex) {
        int[] states = new int[patternCount];

        for (int i = 0; i < patternCount; i++) {

            states[i] = patternList.getPatternState(sequenceIndex, i);
        }

        if (this.states[nodeIndex] == null) {
            this.states[nodeIndex] = new int[patternCount];
        }
        System.arraycopy(states, 0, this.states[nodeIndex], 0, patternCount);

    }


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        updatePartials = true;
        fireModelChanged();
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected void storeState() {
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void restoreState() {
        updatePartials = true;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
    }

    public double[] getTipPartials(int nodeIndex) {
        if (updatePartials) {
            double base = baseDamageRateParameter.getParameterValue(0);
            double factor = ageRateFactorParameter.getParameterValue(0);

            int[] states = this.states[nodeIndex];
            double[] partials = this.partials[nodeIndex];
            for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
                double age = treeModel.getNodeHeight(treeModel.getExternalNode(i));

                double pDamage = base * Math.exp(-factor * age);

                int k = 0;
                for (int j = 0; j < patternCount; j++) {
                    switch (states[j]) {
                        case 0: // is an A
                            partials[k] = 1.0 - pDamage;
                            partials[k + 1] = 0.0;
                            partials[k + 2] = pDamage;
                            partials[k + 3] = 0.0;
                            break;
                        case 1: // is an C
                            partials[k] = 0.0;
                            partials[k + 1] = 1.0 - pDamage;
                            partials[k + 2] = 0.0;
                            partials[k + 3] = pDamage;
                            break;
                        case 2: // is an G
                            partials[k] = pDamage;
                            partials[k + 1] = 0.0;
                            partials[k + 2] = 1.0 - pDamage;
                            partials[k + 3] = 0.0;
                            break;
                        case 3: // is an T
                            partials[k] = 0.0;
                            partials[k + 1] = pDamage;
                            partials[k + 2] = 0.0;
                            partials[k + 3] = 1.0 - pDamage;
                            break;
	                    default: // is an ambiguity
                            partials[k] = 1.0;
                            partials[k + 1] = 1.0;
                            partials[k + 2] = 1.0;
                            partials[k + 3] = 1.0;
                    }
                    k += stateCount;
                }
            }
            updatePartials = false;
        }

        return partials[nodeIndex];
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return ADNA_DAMAGE_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);

            PatternList patternList = (PatternList)xo.getChild(PatternList.class);

            Parameter baseDamageRateParameter = (Parameter)xo.getSocketChild(BASE_DAMAGE_RATE);
            Parameter ageRateFactorParameter = (Parameter)xo.getSocketChild(AGE_RATE_FACTOR);

            ADNADamageModel aDNADamageModel =  new ADNADamageModel(treeModel, patternList, baseDamageRateParameter, ageRateFactorParameter);

            System.out.println("Using aDNA damage model.");

            return aDNADamageModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns a model that allows for post-mortem DNA damage.";
        }

        public Class getReturnType() { return ADNADamageModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),
                new ElementRule(PatternList.class),
                new ElementRule(BASE_DAMAGE_RATE, Parameter.class, "The base rate of accumulation of post-mortem damage", false),
                new ElementRule(AGE_RATE_FACTOR, Parameter.class, "The factor by which rate of damage scales with age of sample", false)
        };
    };

    private final int patternCount;
    private final int stateCount;

    protected int[][] states;
    protected double[][] partials;

    private final TreeModel treeModel;
    private final Parameter baseDamageRateParameter;
    private final Parameter ageRateFactorParameter;

    private boolean updatePartials = true;
}