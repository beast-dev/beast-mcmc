package dr.app.beauti.components;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.PriorType;
import dr.app.beauti.enumTypes.TipDateSamplingType;
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
public class TipDateSamplingComponentOptions implements ComponentOptions {

    private final BeautiOptions options;

    public TipDateSamplingComponentOptions(final BeautiOptions options) {
        this.options = options;
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter("treeModel.tipDates", "date of specified tips", PriorScaleType.TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        modelOptions.createScaleOperator("treeModel.tipDates", modelOptions.demoTuning, 3.0);
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
                    parameter = new Parameter.Builder("age(" + taxon.getId() + ")", "sampled age of taxon, " + taxon.getId())
                            .scaleType(PriorScaleType.TIME_SCALE).prior(PriorType.UNIFORM_PRIOR).initial(height).lower(0.0)
                            .upper(Double.POSITIVE_INFINITY).build();
                    parameter.setPriorEdited(true);
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

    public TaxonList getTaxonSet() {
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