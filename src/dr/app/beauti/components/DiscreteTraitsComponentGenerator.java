package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.enumTypes.TipDateSamplingType;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evoxml.TaxonParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentGenerator extends BaseComponentGenerator {

    public DiscreteTraitsComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        DiscreteTraitsComponentOptions comp = (DiscreteTraitsComponentOptions)options.getComponentOptions(DiscreteTraitsComponentOptions.class);

        if (comp.discreteTraitName == null) {
            return false;
        }

        switch (point) {
            case AFTER_PATTERNS:
            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_PRIOR:
            case IN_MCMC_LIKELIHOOD:
            case IN_FILE_LOG_PARAMETERS:
            case IN_FILE_LOG_LIKELIHOODS:
            case AFTER_FILE_LOG:
            case IN_TREES_LOG:
                return true;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions)options.getComponentOptions(TipDateSamplingComponentOptions.class);

        switch (point) {
            case AFTER_PATTERNS:

            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_PRIOR:
            case IN_MCMC_LIKELIHOOD:
            case IN_FILE_LOG_PARAMETERS:
            case IN_FILE_LOG_LIKELIHOODS:
            case AFTER_FILE_LOG:
            case IN_TREES_LOG:
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Discrete Traits Model";
    }

}