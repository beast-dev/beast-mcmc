package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.substmodel.GeneralSubstitutionModelParser;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AttributePatternsParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.TaxaParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Walter Xie
 */
public class DiscreteTraitGenerator extends Generator {

    public static final String DATA = "DataType";


    public DiscreteTraitGenerator(BeautiOptions options, ComponentFactory[] components) {
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
    public void writeTaxonTraits(Taxon taxon, XMLWriter writer) throws Arguments.ArgumentException {
        for (TraitData trait : options.traits) {
            // may as well write the species traits - may be useful...
//            if (!trait.getName().equalsIgnoreCase(TraitData.TRAIT_SPECIES)) {

//                if (!taxon.containsAttribute(trait.getName())) {
//                    throw new Arguments.ArgumentException("Cannot find trait '" + trait.getName()
//                            + "' in taxon '" + taxon.getId() + "'");
//                }

                writer.writeOpenTag(AttributeParser.ATTRIBUTE, new Attribute[]{
                        new Attribute.Default<String>(Attribute.NAME, trait.getName())});

                // denotes missing data using '?'
                writer.writeText(taxon.containsAttribute(trait.getName()) ? taxon.getAttribute(trait.getName()).toString() : "?");
                writer.writeCloseTag(AttributeParser.ATTRIBUTE);
//            }
        }
    }

    /**
     * write <generalDataType> and <attributePatterns>
     *
     * @param traitData TraitData
     * @param writer    XMLWriter
     */
    public void writeGeneralDataType(TraitData traitData, XMLWriter writer) {
        writer.writeComment("trait = " + traitData.getName() + " trait_type = " + traitData.getTraitType());

        Set<String> states = new HashSet<String>();
        for (PartitionData partition : options.getAllPartitionData(traitData)) {
            states.addAll(partition.getTrait().getStatesOfTrait(options.taxonList));
        }

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitData.getName() + DATA)});

        int numOfStates = states.size();
        writer.writeComment("Number Of States = " + numOfStates);

        for (String eachGD : states) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitData.getName() + AttributePatternsParser.ATTRIBUTE_PATTERNS),
                new Attribute.Default<String>(AttributePatternsParser.ATTRIBUTE, traitData.getName())});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitData.getName() + DATA);
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }


    /**
     * Ancestral Tree Likelihood
     *
     * @param partition PartitionData
     * @param writer    XMLWriter
     */
    public void writeAncestralTreeLikelihood(PartitionData partition, XMLWriter writer) {
        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        writer.writeOpenTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, partition.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD)});

        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, partition.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS);
        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        writer.writeIDref(SiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
        writer.writeIDref(GeneralSubstitutionModelParser.GENERAL_SUBSTITUTION_MODEL, substModel.getPrefix() + AbstractSubstitutionModel.MODEL);

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case UNCORRELATED:
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case RANDOM_LOCAL_CLOCK:
            	writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case AUTOCORRELATED:
                writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        writer.writeCloseTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD);
    }

    public void writeAncestralTreeLikelihoodReferences(XMLWriter writer) {
        for (PartitionData partition : options.dataPartitions) {
            TraitData trait = partition.getTrait();

            if (trait != null && trait.getTraitType() == TraitData.TraitType.DISCRETE) {
                writer.writeIDref(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD,
                        partition.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD);
            }
        }
    }

}

