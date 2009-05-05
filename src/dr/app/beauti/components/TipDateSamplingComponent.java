package dr.app.beauti.components;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TipDateSamplingComponent extends BaseComponentGenerator implements ComponentOptions {

    public TipDateSamplingComponent(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {

        if (tipDateSamplingType == TipDateSamplingType.NO_SAMPLING) {
            return false;
        }

        switch (point) {
            case IN_TREE_MODEL:
            case IN_FILE_LOG_PARAMETERS:
                return true;
            case IN_MCMC_PRIOR:
                return tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final XMLWriter writer) {
        switch (point) {
            case IN_TREE_MODEL: {
                if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                    TaxonList taxa = getTaxonSet();
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeOpenTag("leafHeight",
                                new Attribute[]{
                                        new Attribute.Default<String>("taxon", taxon.getId()),
                                }
                        );
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "age(" + taxon.getId() + ")"), true);
                        writer.writeCloseTag("leafHeight");
                    }
                } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
//                    writer.writeOpenTag("nodeHeights",
//                            new Attribute[]{
//                                    new Attribute.Default<Boolean>("external", taxon.getId()),
//                            }
//                    );
//                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "age(" + taxon.getId() + ")"), true);
//                    writer.writeCloseTag("leafHeight");
                }
            } break;
            case IN_MCMC_PRIOR:
                if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                    // nothing to do - indivual parameter priors are written automatically
                } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {

                }
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                    TaxonList taxa = getTaxonSet();
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "age(" + taxon.getId() + ")"), true);
                    }
                } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, "treeModel.tipDates"), true);
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Tip date sampling";
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter("treeModel.tipDates", "date of specified tips", ModelOptions.TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        modelOptions.createScaleOperator("treeModel.tipDates", 3.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
            TaxonList taxa = getTaxonSet();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                Taxon taxon = taxa.getTaxon(i);
                Parameter parameter = tipDateParameters.get(taxon);
                double height = 0.0;
                if (taxon.getAttribute("height") != null) {
                    height = (Double)taxon.getAttribute("height");
                }
                if (parameter == null) {
                    parameter = new Parameter("age(" + taxon.getId() + ")", 
                            "sampled age of taxon, " + taxon.getId(),
                            ModelOptions.TIME_SCALE,
                            height,
                            0.0, Double.POSITIVE_INFINITY);
                    parameter.priorEdited = true;
                    tipDateParameters.put(taxon, parameter);
                }
                params.add(parameter);
            }
        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
            params.add(modelOptions.getParameter("treeModel.tipDates"));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
            TaxonList taxa = getTaxonSet();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                Taxon taxon = taxa.getTaxon(i);
                Operator operator = tipDateOperators.get(taxon);
                if (operator == null) {
                    Parameter parameter = tipDateParameters.get(taxon);
                    operator = new Operator("age(" + taxon.getId() + ")", "", parameter, OperatorType.SCALE, 0.75, 1.0);                    
                    tipDateOperators.put(taxon, operator);
                }
                ops.add(operator);
            }
        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
            ops.add(modelOptions.getOperator("treeModel.tipDates"));
        }
    }

    private TaxonList getTaxonSet() {
        TaxonList taxa = options.taxonList;

        if (tipDateSamplingTaxonSet != null) {
            taxa = tipDateSamplingTaxonSet;
        }
        return taxa;
    }

    public TipDateSamplingType tipDateSamplingType = TipDateSamplingType.NO_SAMPLING;
    public TaxonList tipDateSamplingTaxonSet = null;

    private Map<Taxon, Parameter> tipDateParameters = new HashMap<Taxon, Parameter>();
    private Map<Taxon, Operator> tipDateOperators = new HashMap<Taxon, Operator>();
}