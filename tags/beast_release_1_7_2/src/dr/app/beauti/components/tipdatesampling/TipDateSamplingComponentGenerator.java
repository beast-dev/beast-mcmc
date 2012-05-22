package dr.app.beauti.components.tipdatesampling;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.types.TipDateSamplingType;
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
public class TipDateSamplingComponentGenerator extends BaseComponentGenerator {

    public TipDateSamplingComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions)options.getComponentOptions(TipDateSamplingComponentOptions.class);

        if (comp.tipDateSamplingType == TipDateSamplingType.NO_SAMPLING) {
            return false;
        }

        switch (point) {
            case IN_TREE_MODEL:
            case IN_FILE_LOG_PARAMETERS:
                return true;
            case AFTER_TREE_MODEL:
            case IN_MCMC_PRIOR:
                return comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions)options.getComponentOptions(TipDateSamplingComponentOptions.class);

        TaxonList taxa = comp.getTaxonSet();
        if (taxa == null) { // all taxa...
            taxa = options.taxonList;
        }

        switch (point) {
            case IN_TREE_MODEL: {
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    Taxon taxon = taxa.getTaxon(i);
                    writer.writeOpenTag("leafHeight",
                            new Attribute[]{
                                    new Attribute.Default<String>(TaxonParser.TAXON, taxon.getId()),
                            }
                    );
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "age(" + taxon.getId() + ")"), true);
                    writer.writeCloseTag("leafHeight");
                }

            } break;
            case AFTER_TREE_MODEL:
                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
                    writer.writeOpenTag("compoundParameter",
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "treeModel.tipDates"),
                            }
                    );
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeIDref(ParameterParser.PARAMETER, "age(" + taxon.getId() + ")");
                    }

                    writer.writeCloseTag("compoundParameter");
                }
                break;
            case IN_MCMC_PRIOR:
                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                    // nothing to do - individual parameter priors are written automatically
                } else if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {

                }
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeIDref(ParameterParser.PARAMETER, "age(" + taxon.getId() + ")");
                    }
                } else if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
                    writer.writeIDref(ParameterParser.PARAMETER, "treeModel.tipDates");
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Tip date sampling";
    }

}