package dr.app.beauti.components.hpm;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.XMLWriter;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public class HierarchicalModelComponentGenerator extends BaseComponentGenerator {

    public HierarchicalModelComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
//        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
//                options.getComponentOptions(HierarchicalModelComponentOptions.class);

        switch (point) {
            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_PRIOR:
                return true;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
//        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
//                options.getComponentOptions(HierarchicalModelComponentOptions.class);

        switch (point) {
            case AFTER_TREE_LIKELIHOOD:
            case IN_MCMC_PRIOR:
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }
    }

    protected String getCommentLabel() {
        return "Hierarchical phylogenetic models";
    }

}