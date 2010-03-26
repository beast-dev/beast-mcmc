package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.DiscreteTraitOptions;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.options.TraitOptions;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.util.Taxon;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.evomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AttributePatternsParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.TaxaParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.model.ProductStatisticParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Walter Xie
 */
public class GeneralTraitGenerator extends Generator {

    public static final String DATA = "Data";


    private int numOfSates; // used in private String getIndicatorsParaValue()

    public GeneralTraitGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * write attr in taxon.
     *
     * @param taxon  Taxon
     * @param writer XMLWriter
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    public void writeAtrrTrait(Taxon taxon, XMLWriter writer) throws Arguments.ArgumentException {
        for (TraitData trait : BeautiOptions.getDiscreteIntegerTraits()) {
            if (!trait.getTraitName().equalsIgnoreCase(TraitOptions.Traits.TRAIT_SPECIES.toString())) {

                if (taxon.containsAttribute(trait.getTraitName())) {
                    throw new Arguments.ArgumentException("Cannot find trait " + trait.getTraitName()
                            + "\nin taxon " + taxon.getId());
                }

                writer.writeOpenTag(Attribute.ATTRIBUTE, new Attribute[]{
                        new Attribute.Default<String>(Attribute.NAME, trait.getTraitName())});

                taxon.getAttribute(trait.getTraitName());
                writer.writeCloseTag(SpeciesBindingsParser.GENE_TREES);
            }
        }
    }

    /**
     * write <generalDataType> and <attributePatterns>
     *
     * @param traitOptions DiscreteTraitOptions
     * @param writer       XMLWriter
     */
    public void writeGeneralDataType(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        List<String> generalData = TraitOptions.getStatesListOfTrait(options.taxonList, traitOptions.getName());

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + DATA)});

        numOfSates = generalData.size();
        writer.writeComment("Number Of Sates = " + numOfSates);

        for (String eachGD : generalData) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + AttributePatternsParser.PATTERNS)});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitOptions.getPrefix() + DATA);
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }

    /**
     * Location Subst Model, Site model, AncestralTreeLikelihood for Discrete Trait
     *
     * @param traitOptions DiscreteTraitOptions
     * @param writer       XMLWriter
     */
    public void writeLocationSubstSiteModel(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        if (traitOptions.getLocationSubstType() == DiscreteTraitOptions.LocationSubstModelType.SYM_SUBST) {
            writer.writeComment("symmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + AbstractSubstitutionModel.MODEL)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitOptions.getName() + DATA);

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeFrequencyModel(traitOptions, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(traitOptions, numOfSates * (numOfSates - 1) / 2, 1, writer);//TODO alway 1?

            writer.writeCloseTag(SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL);

        } else if (traitOptions.getLocationSubstType() == DiscreteTraitOptions.LocationSubstModelType.ASYM_SUBST) {
            writer.writeComment("asymmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + AbstractSubstitutionModel.MODEL),
                    new Attribute.Default<Boolean>(ComplexSubstitutionModelParser.RANDOMIZE, false)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitOptions.getName() + DATA);

            writer.writeOpenTag(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);

            writeFrequencyModel(traitOptions, null, writer);

            writer.writeCloseTag(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(traitOptions, numOfSates * (numOfSates - 1), null, writer);

            writer.writeCloseTag(ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL);

        } else {

        }

        if (traitOptions.isActivateBSSVS()) writeStatisticModel(traitOptions, writer);

        writeSiteModel(traitOptions, writer);
        writeAncestralTreeLikelihood(traitOptions, writer);

    }

    public void writeLogs(DiscreteTraitOptions traitOptions, XMLWriter writer) {

    }

    private void writeFrequencyModel(DiscreteTraitOptions traitOptions, Boolean normalize, XMLWriter writer) {
        if (normalize == null) {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL);
        } else {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<Boolean>(FrequencyModelParser.NORMALIZE, normalize)});
        }

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitOptions.getName() + DATA);

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(traitOptions.getPrefix() + "frequencies", numOfSates, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
    }

    private void writeRatesAndIndicators(DiscreteTraitOptions traitOptions, int dimension, Integer relativeTo, XMLWriter writer) {
        writer.writeComment("rates and indicators");

        if (relativeTo == null) {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES);
        } else {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES, new Attribute[]{
                    new Attribute.Default<Integer>(GeneralSubstitutionModelParser.RELATIVE_TO, relativeTo)});
        }
        traitOptions.getParameter("rates").isFixed = true;
        writeParameter(traitOptions.getParameter("rates"), dimension, writer);
        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (traitOptions.isActivateBSSVS()) {
            writer.writeOpenTag(SVSGeneralSubstitutionModel.INDICATOR);
            traitOptions.getParameter("indicators").isFixed = true;
            writeParameter(traitOptions.getParameter("indicators"), dimension, writer);
            writer.writeCloseTag(SVSGeneralSubstitutionModel.INDICATOR);
        }

    }

    private void writeStatisticModel(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, "nonZeroRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
        writer.writeIDref(ParameterParser.PARAMETER, traitOptions.getPrefix() + "indicators");
        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, "actualRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, false)});
        writer.writeIDref(ParameterParser.PARAMETER, traitOptions.getPrefix() + "indicators");
        writer.writeIDref(ParameterParser.PARAMETER, traitOptions.getPrefix() + "rates");
        writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);
    }

    private void writeSiteModel(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        writer.writeOpenTag(SiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + SiteModel.SITE_MODEL)});

        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writeSubstModelRef(traitOptions, writer);
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        writer.writeOpenTag(GammaSiteModelParser.MUTATION_RATE);
        writeParameter(traitOptions.getParameter("mu"), -1, writer);
        writer.writeCloseTag(GammaSiteModelParser.MUTATION_RATE);

        writer.writeCloseTag(SiteModel.SITE_MODEL);
    }


    private void writeAncestralTreeLikelihood(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        writer.writeOpenTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitOptions.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD)});
        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, traitOptions.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS);
        writer.writeIDref(TreeModel.TREE_MODEL, ""); //TODO
        writer.writeIDref(SiteModel.SITE_MODEL, traitOptions.getPrefix() + SiteModel.SITE_MODEL);
        writeSubstModelRef(traitOptions, writer);
        writer.writeCloseTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD);
    }


    private void writeSubstModelRef(DiscreteTraitOptions traitOptions, XMLWriter writer) {
        String substModel = null;
        if (traitOptions.getLocationSubstType() == DiscreteTraitOptions.LocationSubstModelType.SYM_SUBST) {
            substModel = SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL;
        } else if (traitOptions.getLocationSubstType() == DiscreteTraitOptions.LocationSubstModelType.ASYM_SUBST) {
            substModel = ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL;
        } else {

        }
        writer.writeIDref(substModel, traitOptions.getPrefix() + AbstractSubstitutionModel.MODEL);
    }


}

