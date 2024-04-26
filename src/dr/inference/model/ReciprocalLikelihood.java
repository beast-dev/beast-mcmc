package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.xml.Reportable;

import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class ReciprocalLikelihood implements Likelihood, Reportable {

    private final Likelihood likelihood;

    public ReciprocalLikelihood(Likelihood likelihood) {
        this.likelihood = likelihood;
    }

    @Override
    public LogColumn[] getColumns() {
        return likelihood.getColumns();
    }

    @Override
    public Model getModel() {
        return likelihood.getModel();
    }

    @Override
    public double getLogLikelihood() {
        return -likelihood.getLogLikelihood();
    }

    @Override
    public void makeDirty() {
        likelihood.makeDirty();
    }

    @Override
    public String prettyName() {
        return likelihood.prettyName();
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return likelihood.getLikelihoodSet();
    }

    @Override
    public boolean isUsed() {
        return likelihood.isUsed();
    }

    @Override
    public void setUsed() {
        likelihood.setUsed();
    }

    @Override
    public boolean evaluateEarly() {
        return likelihood.evaluateEarly();
    }

    @Override
    public String getId() {
        return likelihood.getId();
    }

    @Override
    public void setId(String id) {
        likelihood.setId(id);
    }

    @Override
    public String getReport() {
        return "Reciprocal likelihood " + getId() + " " + getLogLikelihood();
    }
}
