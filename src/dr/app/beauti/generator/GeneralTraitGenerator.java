package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.LocationSubstModelType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.util.Taxon;
import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.evomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.util.Attribute;

/**
 * @author Walter Xie
 */
public class GeneralTraitGenerator extends Generator {

    public static final String DATA = "Data";


    private int numOfSates; // used in private String getIndicatorsParaValue()

    public GeneralTraitGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    public static String getLocationSubstModelTag(PartitionSubstitutionModel substModel) {
        if (substModel.getLocationSubstType() == LocationSubstModelType.SYM_SUBST) {
            return SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL;
        } else if (substModel.getLocationSubstType() == LocationSubstModelType.ASYM_SUBST) {
            return ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL;
        } else {
            return null;
        }
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
            if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString())) {

                if (taxon.containsAttribute(trait.getName())) {
                    throw new Arguments.ArgumentException("Cannot find trait " + trait.getName()
                            + "\nin taxon " + taxon.getId());
                }

                writer.writeOpenTag(Attribute.ATTRIBUTE, new Attribute[]{
                        new Attribute.Default<String>(Attribute.NAME, trait.getName())});

                taxon.getAttribute(trait.getName());
                writer.writeCloseTag(SpeciesBindingsParser.GENE_TREES);
            }
        }
    }

    /**
     * write <generalDataType> and <attributePatterns>
     *
     * @param TraitData DiscreteTraitData
     * @param writer       XMLWriter
     */
/*    public void writeGeneralDataType(DiscreteTraitData TraitData, XMLWriter writer) {
        List<String> generalData = TraitData.getStatesListOfTrait(BeautiOptions.taxonList, TraitData.getName());

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + DATA)});

        numOfSates = generalData.size();
        writer.writeComment("Number Of Sates = " + numOfSates);

        for (String eachGD : generalData) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + AttributePatternsParser.PATTERNS)});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, TraitData.getPrefix() + DATA);
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }

    /**
     * Location Subst Model, Site model, AncestralTreeLikelihood for Discrete Trait
     *
     * @param TraitData DiscreteTraitData
     * @param writer       XMLWriter
     
    public void writeLocationSubstSiteModel(DiscreteTraitData TraitData, XMLWriter writer) {
        if (TraitData.getLocationSubstType() == LocationSubstModelType.SYM_SUBST) {
            writer.writeComment("symmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + AbstractSubstitutionModel.MODEL)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, TraitData.getName() + DATA);

            writer.writeOpenTag(GeneralSubstitutionModelParser.FREQUENCIES);

            writeFrequencyModel(TraitData, true, writer);

            writer.writeCloseTag(GeneralSubstitutionModelParser.FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(TraitData, numOfSates * (numOfSates - 1) / 2, 1, writer);//TODO alway 1?

            writer.writeCloseTag(SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL);

        } else if (TraitData.getLocationSubstType() == LocationSubstModelType.ASYM_SUBST) {
            writer.writeComment("asymmetric CTMC model for discrete state reconstructions");

            writer.writeOpenTag(ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + AbstractSubstitutionModel.MODEL),
                    new Attribute.Default<Boolean>(ComplexSubstitutionModelParser.RANDOMIZE, false)});

            writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, TraitData.getName() + DATA);

            writer.writeOpenTag(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);

            writeFrequencyModel(TraitData, null, writer);

            writer.writeCloseTag(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);

            //---------------- rates and indicators -----------------
            writeRatesAndIndicators(TraitData, numOfSates * (numOfSates - 1), null, writer);

            writer.writeCloseTag(ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL);

        } else {

        }

        if (TraitData.isActivateBSSVS()) writeStatisticModel(TraitData, writer);

        writeSiteModel(TraitData, writer);
        writeAncestralTreeLikelihood(TraitData, writer);

    }

    public void writeLogs(DiscreteTraitData TraitData, XMLWriter writer) {

    }

    private void writeFrequencyModel(DiscreteTraitData TraitData, Boolean normalize, XMLWriter writer) {
        if (normalize == null) {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL);
        } else {
            writer.writeOpenTag(FrequencyModelParser.FREQUENCY_MODEL, new Attribute[]{
                    new Attribute.Default<Boolean>(FrequencyModelParser.NORMALIZE, normalize)});
        }

        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, TraitData.getName() + DATA);

        writer.writeOpenTag(FrequencyModelParser.FREQUENCIES);
        writeParameter(TraitData.getPrefix() + "frequencies", numOfSates, Double.NaN, Double.NaN, Double.NaN, writer);
        writer.writeCloseTag(FrequencyModelParser.FREQUENCIES);

        writer.writeCloseTag(FrequencyModelParser.FREQUENCY_MODEL);
    }

    private void writeRatesAndIndicators(DiscreteTraitData TraitData, int dimension, Integer relativeTo, XMLWriter writer) {
        writer.writeComment("rates and indicators");

        if (relativeTo == null) {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES);
        } else {
            writer.writeOpenTag(GeneralSubstitutionModelParser.RATES, new Attribute[]{
                    new Attribute.Default<Integer>(GeneralSubstitutionModelParser.RELATIVE_TO, relativeTo)});
        }
        TraitData.getParameter("rates").isFixed = true;
        writeParameter(TraitData.getParameter("rates"), dimension, writer);
        writer.writeCloseTag(GeneralSubstitutionModelParser.RATES);

        if (TraitData.isActivateBSSVS()) {
            writer.writeOpenTag(SVSGeneralSubstitutionModel.INDICATOR);
            TraitData.getParameter("indicators").isFixed = true;
            writeParameter(TraitData.getParameter("indicators"), dimension, writer);
            writer.writeCloseTag(SVSGeneralSubstitutionModel.INDICATOR);
        }

    }

    private void writeStatisticModel(DiscreteTraitData TraitData, XMLWriter writer) {
        writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, "nonZeroRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
        writer.writeIDref(ParameterParser.PARAMETER, TraitData.getPrefix() + "indicators");
        writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);

        writer.writeOpenTag(ProductStatisticParser.PRODUCT_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, "actualRates"),
                new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, false)});
        writer.writeIDref(ParameterParser.PARAMETER, TraitData.getPrefix() + "indicators");
        writer.writeIDref(ParameterParser.PARAMETER, TraitData.getPrefix() + "rates");
        writer.writeCloseTag(ProductStatisticParser.PRODUCT_STATISTIC);
    }

    private void writeSiteModel(DiscreteTraitData TraitData, XMLWriter writer) {
        writer.writeOpenTag(SiteModel.SITE_MODEL, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + SiteModel.SITE_MODEL)});

        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writeSubstModelRef(TraitData, writer);
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        writer.writeOpenTag(GammaSiteModelParser.MUTATION_RATE);
        writeParameter(TraitData.getParameter("mu"), -1, writer);
        writer.writeCloseTag(GammaSiteModelParser.MUTATION_RATE);

        writer.writeCloseTag(SiteModel.SITE_MODEL);
    }


    private void writeAncestralTreeLikelihood(DiscreteTraitData TraitData, XMLWriter writer) {
        writer.writeOpenTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TraitData.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD)});
        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, TraitData.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS);
        writer.writeIDref(TreeModel.TREE_MODEL, ""); //TODO
        writer.writeIDref(SiteModel.SITE_MODEL, TraitData.getPrefix() + SiteModel.SITE_MODEL);
        writeSubstModelRef(TraitData, writer);
        writer.writeCloseTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD);
    }

*/



}

