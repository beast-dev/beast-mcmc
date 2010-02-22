package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.MutationDeathModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Package: ALSTreeLikelihood
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 13, 2008
 * Time: 10:13:07 AM
 */
public class ALSTreeLikelihood extends TreeLikelihood {
    public static final String LIKE_NAME = "alsTreeLikelihood";
    public static final String INTEGRATE_GAIN_RATE = "integrateGainRate";
    public static final String OBSERVATION_PROCESS = "observationProcess";
    public static final String OBSERVATION_TYPE = "type";
    public static final String OBSERVATION_TAXON = "taxon";
    final static String IMMIGRATION_RATE = "immigrationRate";

    protected AbstractObservationProcess observationProcess;

    public ALSTreeLikelihood(AbstractObservationProcess observationProcess, PatternList patternList, TreeModel treeModel, SiteModel siteModel, BranchRateModel branchRateModel, boolean useAmbiguities, boolean storePartials) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials, false, false);

        this.observationProcess = observationProcess;
        addModel(observationProcess);

        // TreeLikelihood does not initialize the partials for tips, we'll do it ourselves
        int extNodeCount = treeModel.getExternalNodeCount();
        for (int i = 0; i < extNodeCount; i++) {
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            setPartials(likelihoodCore, patternList, categoryCount, index, i);
        }

/*        //Examine the tree
        double totalTime=0.0;
        double realTime = 0.0;
        for(int i=0; i<treeModel.getNodeCount();++i){
            NodeRef node = treeModel.getNode(i);
            double branchRate = branchRateModel.getBranchRate(treeModel,node);
            double branchTime = treeModel.getBranchLength(node);
            totalTime+=branchRate*branchTime;
            realTime += branchTime;
            System.err.println("Node "+node.toString()+ " time: "+branchTime+ " rate "+branchRate+" together "+branchTime*branchRate);
        }
        System.err.println("TotalTime: "+totalTime);
        System.err.println("RealTime: "+realTime);*/

    }

    protected double calculateLogLikelihood() {
        // Calculate the partial likelihoods
        super.calculateLogLikelihood();
        // get the frequency model
        double[] freqs = frequencyModel.getFrequencies();
        // let the observationProcess handle the rest
        return observationProcess.nodePatternLikelihood(freqs, likelihoodCore);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == observationProcess || model == treeModel) {
            makeDirty();
        } else
            super.handleModelChangedEvent(model, object, index);
    }

    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LIKE_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = false;
            boolean storePartials = true;
            if (xo.hasAttribute(TreeLikelihoodParser.USE_AMBIGUITIES)) {
                useAmbiguities = xo.getBooleanAttribute(TreeLikelihoodParser.USE_AMBIGUITIES);
            }
            if (xo.hasAttribute(TreeLikelihoodParser.STORE_PARTIALS)) {
                storePartials = xo.getBooleanAttribute(TreeLikelihoodParser.STORE_PARTIALS);
            }

            boolean integrateGainRate = xo.getBooleanAttribute(INTEGRATE_GAIN_RATE);

            //AbstractObservationProcess observationProcess = (AbstractObservationProcess) xo.getChild(AbstractObservationProcess.class);


            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
            Parameter mu = ((MutationDeathModel) siteModel.getSubstitutionModel()).getDeathParameter();
            Parameter lam;
            if (!integrateGainRate) {
                lam = (Parameter) xo.getElementFirstChild(IMMIGRATION_RATE);
            } else {
                lam = new Parameter.Default("gainRate", 1.0, 0.001, 1.999);
            }
            AbstractObservationProcess observationProcess = null;

            Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating ALSTreeLikelihood model.");
            for (int i = 0; i < xo.getChildCount(); ++i) {
                Object cxo = xo.getChild(i);
                if (cxo instanceof XMLObject && ((XMLObject) cxo).getName().equals(OBSERVATION_PROCESS)) {
                    if (((XMLObject) cxo).getStringAttribute(OBSERVATION_TYPE).equals("singleTip")) {
                        String taxonName = ((XMLObject) cxo).getStringAttribute(OBSERVATION_TAXON);
                        Taxon taxon = treeModel.getTaxon(treeModel.getTaxonIndex(taxonName));
                        observationProcess = new SingleTipObservationProcess(treeModel, patternList, siteModel,
                                branchRateModel, mu, lam, taxon);
                        Logger.getLogger("dr.evolution").info("All traits are assumed extant in " + taxonName);
                    } else {  // "anyTip" observation process
                        observationProcess = new AnyTipObservationProcess("anyTip", treeModel, patternList,
                                siteModel, branchRateModel, mu, lam);
                        Logger.getLogger("dr.evolution").info("Observed traits are assumed to be extant in at least one tip node.");
                    }

                    observationProcess.setIntegrateGainRate(integrateGainRate);
                }
            }
            Logger.getLogger("dr.evolution").info("\tIf you publish results using Acquisition-Loss-Mutaion (ALS) Model likelihood, please reference Alekseyenko, Lee and Suchard (2008) Syst. Biol 57: 772-784.\n---------------------------------\n");

            return new ALSTreeLikelihood(observationProcess, patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of a patternlist on a tree given the site model.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
                AttributeRule.newBooleanRule(TreeLikelihoodParser.STORE_PARTIALS, true),
                AttributeRule.newBooleanRule(INTEGRATE_GAIN_RATE),
                new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(OBSERVATION_PROCESS,
                        new XMLSyntaxRule[]{AttributeRule.newStringRule(OBSERVATION_TYPE, false),
                                AttributeRule.newStringRule(OBSERVATION_TAXON, true)})
                //new ElementRule(AbstractObservationProcess.class)
        };
    };

}
