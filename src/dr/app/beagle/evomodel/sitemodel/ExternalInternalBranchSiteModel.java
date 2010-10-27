package dr.app.beagle.evomodel.sitemodel;

import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public class ExternalInternalBranchSiteModel extends AbstractModel implements BranchSiteModel, Citable {
    public ExternalInternalBranchSiteModel(List<SubstitutionModel> substModelList, List<FrequencyModel> frequencyModelList) {
        super("ExternalInternalBranchSiteModel");

        if (substModelList.size() != 2) {
            throw new IllegalArgumentException("ExternalInternalBranchSiteModel requires two SubstitutionModels");
        }

        if (frequencyModelList.size() != 1) {
            throw new IllegalArgumentException("ExternalInternalBranchSiteModel requires one FrequencyModel");
        }

        this.substModelList = substModelList;
        this.frequencyModelList = frequencyModelList;

        for (SubstitutionModel model : substModelList) {
            addModel(model);
        }
        for (FrequencyModel model : frequencyModelList) {
            addModel(model);
        }
    }

    public int getBranchIndex(final Tree tree, final NodeRef node) {
        return (tree.isExternal(node) ? 1 : 0);
    }

    public EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex) {
        return substModelList.get(branchIndex).getEigenDecomposition();
    }

    public double[] getStateFrequencies(int categoryIndex) {
        return frequencyModelList.get(categoryIndex).getFrequencies();
    }

    public boolean canReturnComplexDiagonalization() {
        for (SubstitutionModel model : substModelList) {
            if (model.canReturnComplexDiagonalization()) {
                return true;
            }
        }
        return false;
    }
        
    public int getEigenCount() {
        return 2;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    private final List<SubstitutionModel> substModelList;
    private final List<FrequencyModel> frequencyModelList;

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("P", "Lemey"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
        return citations;
    }
}