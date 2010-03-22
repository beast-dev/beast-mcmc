package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.options.TraitsOptions;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.util.Taxon;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.evoxml.AttributePatternsParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.TaxaParser;
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
     * @param taxon
     * @param writer XMLWriter
     */
    public void writeAtrrTrait(Taxon taxon, XMLWriter writer) throws Arguments.ArgumentException {
        for (TraitGuesser trait : TraitsOptions.traits) {
            if (!trait.getTraitName().equalsIgnoreCase(TraitsOptions.Traits.TRAIT_SPECIES.toString())) {

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
     * @param traitName
     * @param writer    XMLWriter
     */
    public void writeGeneralDataType(String traitName, XMLWriter writer) {
        List<String> generalData = TraitsOptions.getStatesListOfTrait(options.taxonList, traitName);

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitName + DATA)});

        numOfSates = generalData.size();
        writer.writeComment("Number Of Sates = " + numOfSates);

        for (String eachGD : generalData) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitName + AttributePatternsParser.PATTERNS)});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitName + DATA);
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }
    
    public void writeLocationSubstModel(XMLWriter writer) {
        writer.writeComment("Collection of Gene Trees");

        writer.writeOpenTag(SpeciesBindingsParser.GENE_TREES, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SpeciesBindingsParser.GENE_TREES)});


        writer.writeCloseTag(SpeciesBindingsParser.GENE_TREES);
    }

    public void writeGeneralSubstModel(XMLWriter writer) {
        writer.writeComment("");


    }

    private void writeSiteModel(XMLWriter writer) {
        writer.writeComment("");



    }


    private void writeAncestralTreeLikelihood(XMLWriter writer) {
        writer.writeComment("");

    }


    private String getIndicatorsParaValue() {
        String v = "";

        return v;
    }
}

