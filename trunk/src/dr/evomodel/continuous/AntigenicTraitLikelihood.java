package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicTraitLikelihood extends AbstractModelLikelihood {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    public AntigenicTraitLikelihood(
            String traitName,
            TreeModel tree,
            MultivariateDiffusionModel diffusionModel,
            CompoundParameter tipTraitParameter,
            CompoundParameter virusTraitParameter,
            CompoundParameter serumTraitParameter,
            double[][] assayTable,
            String[] virusNames,
            String[] serumNames) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        this.virusNames = virusNames;
        this.serumNames = serumNames;

        int tipCount = tipTraitParameter.getDimension();

        // the total number of viruses is the number of rows in the table
        int virusCount = assayTable.length;
        // the number of sera is the number of columns
        int serumCount = assayTable[0].length;

        if (assayTable.length != virusNames.length) {
            throw new IllegalArgumentException("The number of rows in the assay table doesn't match the virus name list");
        }

        if (assayTable[0].length != serumNames.length) {
            throw new IllegalArgumentException("The number of columns in the assay table doesn't match the serum name list");
        }

        this.assayTable = new double[virusCount][serumCount];

        // the tip -> virus map
        tipIndices = new int[tree.getExternalNodeCount()];

        Map<String, Integer> tipNameMap = new HashMap<String, Integer>();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            tipNameMap.put(taxon.getId(), i);

            tipIndices[i] = -1;
        }

        // the virus -> tip map
        virusIndices = new int[virusCount];

        for (int i = 0; i < virusCount; i++) {
            virusIndices[i] = -1;

            // if the virus is in the tree then add a entry to map tip to virus
            Integer tipIndex = tipNameMap.get(virusNames[i]);
            if (tipIndex != null) {
                tipIndices[tipIndex] = i;
                virusIndices[i] = tipIndex;
            } else {
                System.err.println("Virus, " + virusNames[i] + ", not found in tree");
            }

            for (int j = 0; j < serumCount; i++) {
                this.assayTable[i][j] = transform(assayTable[i][j]);
            }
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            if (tipIndices[i] == -1) {
                System.err.println("Tree tip, " + taxon.getId() + ", not found in virus assay table");
            }
        }

        virusTraitParameter.setDimension(virusCount);
        addVariable(virusTraitParameter);

        serumTraitParameter.setDimension(serumCount);
        addVariable(serumTraitParameter);

        // we don't listen to the tip trait parameter because we are setting that
    }

    private double transform(final double value) {
        // transform to log_2
        return Math.log(value) / Math.log(2.0);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return 0;
    }

    public void makeDirty() {
    }


    private final double[][] assayTable;
    private final String[] virusNames;
    private final String[] serumNames;

    private final int[] tipIndices;
    private final int[] virusIndices;
}
