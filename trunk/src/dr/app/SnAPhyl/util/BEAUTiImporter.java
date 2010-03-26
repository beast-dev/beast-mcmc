package dr.app.SnAPhyl.util;

import dr.app.beauti.options.*;
import dr.app.beauti.util.NexusApplicationImporter;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.NexusImporter.NexusBlock;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import org.jdom.JDOMException;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;



/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class BEAUTiImporter {

    private final BeautiOptions options;

    public BEAUTiImporter(BeautiOptions options) {
        this.options = options;
    }
 
    public void importFromFile(File file) throws Exception {
        try {
            Reader reader = new FileReader(file);

            BufferedReader bufferedReader = new BufferedReader(reader);
            String line = bufferedReader.readLine();
            while (line != null && line.length() == 0) {
                line = bufferedReader.readLine();
            }

            if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
                // is a NEXUS file
                importNexusFile(file);
            } else {
                // assume it is a BEAST XML file and see if that works...
                importBEASTFile(file);
            }
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }
    
    // xml
    private void importBEASTFile(File file) throws IOException, JDOMException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            BeastImporter importer = new BeastImporter(reader);

            List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            List<Alignment> alignments = new ArrayList<Alignment>();

            importer.importBEAST(taxonLists, alignments);

            TaxonList taxa = taxonLists.get(0);

            for (Alignment alignment : alignments) {
                setData(taxa, alignment, null, null, null, file.getName());
            }
        } catch (JDOMException e) {
            throw new JDOMException (e.getMessage());
        } catch (ImportException e) {
            throw new ImportException (e.getMessage());
        } catch (IOException e) {
            throw new IOException (e.getMessage());
        }

    }
    
    // nexus
    private void importNexusFile(File file) throws Exception {
        TaxonList taxa = null;
        SimpleAlignment alignment = null;
        List<Tree> trees = new ArrayList<Tree>();
        PartitionSubstitutionModel model = null;
        List<NexusApplicationImporter.CharSet> charSets = new ArrayList<NexusApplicationImporter.CharSet>();

        try {
            FileReader reader = new FileReader(file);


            NexusApplicationImporter importer = new NexusApplicationImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxa != null) {
                            throw new MissingBlockException("TAXA block already defined");
                        }

                        taxa = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.CALIBRATION_BLOCK) {
                        if (taxa == null) {
                            throw new MissingBlockException("TAXA or DATA block must be defined before a CALIBRATION block");
                        }

                        importer.parseCalibrationBlock(taxa);

                    } else if (block == NexusImporter.CHARACTERS_BLOCK) {

                        if (taxa == null) {
                            throw new MissingBlockException("TAXA block must be defined before a CHARACTERS block");
                        }

                        if (alignment != null) {
                            throw new MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        alignment = (SimpleAlignment) importer.parseCharactersBlock(taxa);

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        if (alignment != null) {
                            throw new MissingBlockException("CHARACTERS or DATA block already defined");
                        }

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = (SimpleAlignment) importer.parseDataBlock(taxa);
                        if (taxa == null) {
                            taxa = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        // I guess there is no reason not to allow multiple trees blocks
//                        if (trees.size() > 0) {
//                            throw new MissingBlockException("TREES block already defined");
//                        }

                        Tree[] treeArray = importer.parseTreesBlock(taxa);
                        trees.addAll(Arrays.asList(treeArray));

                        if (taxa == null && trees.size() > 0) {
                            taxa = trees.get(0);
                        }


                    } else if (block == NexusApplicationImporter.PAUP_BLOCK) {

                        model = importer.parsePAUPBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.MRBAYES_BLOCK) {

                        model = importer.parseMrBayesBlock(options, charSets);

                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (IOException e) {
            throw new IOException (e.getMessage());
        } catch (ImportException e) {
            throw new ImportException (e.getMessage());
        } catch (Exception e) {
            throw new Exception (e.getMessage());
        }

        setData(taxa, alignment, trees, model, charSets, file.getName());
    }
    
    //TODO need refactory to simplify
    private void setData(TaxonList taxa, Alignment alignment, List<Tree> trees, PartitionSubstitutionModel model,
            List<NexusApplicationImporter.CharSet> charSets, String fileName) throws ImportException {
        String fileNameStem = dr.app.util.Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "TRE", "TREE", "XML"});

        if (options.taxonList == null) {
            // This is the first partition to be loaded...

            options.taxonList = new Taxa(taxa);

            // check the taxon names for invalid characters
            boolean foundAmp = false;
            for (Taxon taxon : taxa) {
                String name = taxon.getId();
                if (name.indexOf('&') >= 0) {
                    foundAmp = true;
                }
            }
            if (foundAmp) {
                throw new ImportException ("One or more taxon names include an illegal character ('&').\n" +
                        "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                        "Please edit the taxon name(s) before reloading the data file.");
            }

            // make sure they all have dates...
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                if (taxa.getTaxonAttribute(i, "date") == null) {
                    Date origin = new Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxa.getTaxon(i).setAttribute("date", date);
                }
            }

            options.fileNameStem = fileNameStem;

        } else {
            // This is an additional partition so check it uses the same taxa
            if (!options.allowDifferentTaxa) { // not allow Different Taxa
                List<String> oldTaxa = new ArrayList<String>();
                for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                    oldTaxa.add(options.taxonList.getTaxon(i).getId());
                }
                List<String> newTaxa = new ArrayList<String>();
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    newTaxa.add(taxa.getTaxon(i).getId());
                }

                if (!(oldTaxa.containsAll(newTaxa) && oldTaxa.size() == newTaxa.size())) {
                    
                    
                }
            } else { // allow Different Taxa
                // AR - it will be much simpler just to consider options.taxonList
                // to be the union set of all taxa. Each data partition has an alignment
                // which is a taxon list containing the taxa specific to that partition

                // add the new diff taxa
                List<String> prevTaxa = new ArrayList<String>();
                for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
                    prevTaxa.add(options.taxonList.getTaxon(i).getId());
                }
                for (int i = 0; i < taxa.getTaxonCount(); i++) {
                    if (!prevTaxa.contains(taxa.getTaxon(i).getId())) {
                        options.taxonList.addTaxon(taxa.getTaxon(i));
                    }
                }

            }
        }

        addAlignment(alignment, charSets, model, fileName, fileNameStem);

        addTrees(trees);
    }

    private void addAlignment(Alignment alignment, List<NexusApplicationImporter.CharSet> charSets,
                              PartitionSubstitutionModel model,
                              String fileName, String fileNameStem) {
        if (alignment != null) {
            List<PartitionData> partitions = new ArrayList<PartitionData>();
            if (charSets != null && charSets.size() > 0) {
                for (NexusApplicationImporter.CharSet charSet : charSets) {
                    partitions.add(new PartitionData(charSet.getName(), fileName,
                            alignment, charSet.getFromSite(), charSet.getToSite(), charSet.getEvery()));
                }
            } else {
                partitions.add(new PartitionData(fileNameStem, fileName, alignment));
            }
            for (PartitionData partition : partitions) {
                BeautiOptions.dataPartitions.add(partition);
                
                if (model != null) {//TODO Cannot load Clock Model and Tree Model from BEAST file yet                
                    partition.setPartitionSubstitutionModel(model);
                    model.addPartitionData(partition);
                    
                    // use same tree model and same tree prior in beginning
                    if (options.getPartitionTreeModels() != null
                            && options.getPartitionTreeModels().size() == 1) {
                        PartitionTreeModel ptm = options.getPartitionTreeModels().get(0);
                        partition.setPartitionTreeModel(ptm); // same tree model, therefore same prior       
                    }
                    if (partition.getPartitionTreeModel() == null) {
                        // PartitionTreeModel based on PartitionData
                        PartitionTreeModel ptm = new PartitionTreeModel(options, partition);
                        partition.setPartitionTreeModel(ptm);
                        
                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(options, ptm);
                        ptm.setPartitionTreePrior(ptp);
                        
//                        options.addPartitionTreeModel(ptm);
//                        options.shareSameTreePrior = true;                                                
                    }

                    // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                    if (options.getPartitionClockModels() != null
                            && options.getPartitionClockModels().size() == 1) {
                        PartitionClockModel pcm = options.getPartitionClockModels().get(0);                      
                        partition.setPartitionClockModel(pcm);                        
                    }
                    if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                        PartitionClockModel pcm = new PartitionClockModel(options, partition);
                        partition.setPartitionClockModel(pcm);
//                        options.addPartitionClockModel(pcm);
                    }
//                    options.clockModelOptions.fixRateOfFirstClockPartition();

                } else {// only this works                    
                    if (options.getPartitionSubstitutionModels() != null
                            && options.getPartitionSubstitutionModels().size() == 1) { // use same substitution model in beginning
                        PartitionSubstitutionModel psm = options.getPartitionSubstitutionModels().get(0);
                        if (psm.getDataType() == alignment.getDataType()) {
                            partition.setPartitionSubstitutionModel(psm);
                        } else {
                            //TODO exception?
                        }
                    }
                    if (partition.getPartitionSubstitutionModel() == null) {
                        // PartitionSubstitutionModel based on PartitionData
                        PartitionSubstitutionModel psm = new PartitionSubstitutionModel(options, partition);
                        partition.setPartitionSubstitutionModel(psm);
//                        options.addPartitionSubstitutionModel(psm);
                    }

                    // use same tree model and same tree prior in beginning
                    if (options.getPartitionTreeModels() != null
                            && options.getPartitionTreeModels().size() == 1) {
                        PartitionTreeModel ptm = options.getPartitionTreeModels().get(0);
                        partition.setPartitionTreeModel(ptm); // same tree model, therefore same prior
                    }
                    if (partition.getPartitionTreeModel() == null) {
                        // PartitionTreeModel based on PartitionData
                        PartitionTreeModel ptm = new PartitionTreeModel(options, partition);
                        partition.setPartitionTreeModel(ptm);

                        // PartitionTreePrior always based on PartitionTreeModel
                        PartitionTreePrior ptp = new PartitionTreePrior(options, ptm);
                        ptm.setPartitionTreePrior(ptp);

//                        options.addPartitionTreeModel(ptm);
//                        options.shareSameTreePrior = true;                                                
                    }

                    // use same clock model in beginning, have to create after partition.setPartitionTreeModel(ptm);
                    if (options.getPartitionClockModels() != null
                            && options.getPartitionClockModels().size() == 1) {
                        PartitionClockModel pcm = options.getPartitionClockModels().get(0);
                        partition.setPartitionClockModel(pcm);
                    }
                    if (partition.getPartitionClockModel() == null) {
                        // PartitionClockModel based on PartitionData
                        PartitionClockModel pcm = new PartitionClockModel(options, partition);
                        partition.setPartitionClockModel(pcm);
//                        options.addPartitionClockModel(pcm);
                    }
                }
            }
            
            options.updateLinksBetweenPDPCMPSMPTMPTPP();
            options.updatePartitionClockTreeLinks();
            options.clockModelOptions.fixRateOfFirstClockPartition();
        }
    }
    
    private void addTrees(List<Tree> trees) {
        if (trees != null && trees.size() > 0) {
            for (Tree tree : trees) {
                String id = tree.getId();
                if (id == null || id.trim().length() == 0) {
                    tree.setId("tree_" + (options.userTrees.size() + 1));
                } else {
                    String newId = id;
                    int count = 1;
                    for (Tree tree1 : options.userTrees) {
                        if (tree1.getId().equals(newId)) {
                            newId = id + "_" + count;
                            count++;
                        }
                    }
                    tree.setId(newId);
                }
                options.userTrees.add(tree);
            }
        }
    }

}
