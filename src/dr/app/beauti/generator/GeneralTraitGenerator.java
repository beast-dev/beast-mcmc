package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.PopulationSizeModelType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.TraitGuesser;
import dr.app.beauti.options.TraitsOptions;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.evomodelxml.speciation.SpeciesBindingsParser;
import dr.evomodelxml.speciation.SpeciesBindingsSPinfoParser;
import dr.evomodelxml.speciation.SpeciesTreeModelParser;
import dr.evoxml.TaxonParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Walter Xie
 */
public class GeneralTraitGenerator extends Generator {

    private int numOfSpecies; // used in private String getIndicatorsParaValue()

    public GeneralTraitGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * write tag <sp>
     *
     * @param taxonList TaxonList
     * @param writer    XMLWriter
     */
    public void writeMultiSpecies(TaxonList taxonList, XMLWriter writer) {
        List<String> species = options.starBEASTOptions.getSpeciesList();
        String sp;

        numOfSpecies = species.size(); // used in private String getIndicatorsParaValue()

        for (String eachSp : species) {
            writer.writeOpenTag(SpeciesBindingsSPinfoParser.SP, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, eachSp)});

            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                sp = taxon.getAttribute(TraitsOptions.Traits.TRAIT_SPECIES.toString()).toString();

                if (sp.equals(eachSp)) {
                    writer.writeIDref(TaxonParser.TAXON, taxon.getId());
                }

            }
            writer.writeCloseTag(SpeciesBindingsSPinfoParser.SP);
        }

        writeGeneTrees(writer);
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


    public void writeGeneTrees(XMLWriter writer) {
        writer.writeComment("Collection of Gene Trees");

        writer.writeOpenTag(SpeciesBindingsParser.GENE_TREES, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, SpeciesBindingsParser.GENE_TREES)});


        writer.writeCloseTag(SpeciesBindingsParser.GENE_TREES);
    }


    private void writeGeneralDataType(XMLWriter writer) {
        writer.writeComment("Species Tree: Provides Per branch demographic function");

        List<Attribute> attributes = new ArrayList<Attribute>();

        attributes.add(new Attribute.Default<String>(XMLParser.ID, SP_TREE));
        // *BEAST always share same tree prior
        if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS_CONSTANT) {
            attributes.add(new Attribute.Default<String>(SpeciesTreeModelParser.CONST_ROOT_POPULATION, "true"));
        } else if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONSTANT) {
            attributes.add(new Attribute.Default<String>(SpeciesTreeModelParser.CONSTANT_POPULATION, "true"));
        }

        writer.writeOpenTag(SpeciesTreeModelParser.SPECIES_TREE, attributes);

        writer.writeIDref(TraitsOptions.Traits.TRAIT_SPECIES.toString(), TraitsOptions.Traits.TRAIT_SPECIES.toString());

        // take sppSplitPopulations value from partionModel(?).constant.popSize
        // *BEAST always share same tree prior
        double popSizeValue = options.getPartitionTreePriors().get(0).getParameter("constant.popSize").initial; // "initial" is "value"
        writer.writeOpenTag(SpeciesTreeModelParser.SPP_SPLIT_POPULATIONS, new Attribute[]{
                new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(popSizeValue))});

        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SpeciesTreeModelParser.SPECIES_TREE + "." + SPLIT_POPS)}, true);

        writer.writeCloseTag(SpeciesTreeModelParser.SPP_SPLIT_POPULATIONS);

        writer.writeCloseTag(SpeciesTreeModelParser.SPECIES_TREE);

    }

    private void writeSpeciesTreeModel(XMLWriter writer) {
        Parameter para;

        writer.writeComment("Species Tree: tree prior");


        throw new IllegalArgumentException("Get wrong species tree prior using *BEAST : " + options.getPartitionTreePriors().get(0).getNodeHeightPrior().toString());


    }


    private void writeSpeciesTreeLikelihood(XMLWriter writer) {
        writer.writeComment("Species Tree: Likelihood of species tree");


        writer.writeOpenTag(SpeciesTreeModelParser.SPECIES_TREE);
        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);
        writer.writeCloseTag(SpeciesTreeModelParser.SPECIES_TREE);

        writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);
    }


    private String getIndicatorsParaValue() {
        String v = "";

        return v;
    }
}

