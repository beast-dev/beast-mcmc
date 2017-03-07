/*
 * DataModelImporter.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beastgen;

import dr.app.beauti.options.*;
import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.FastaImporter;
import dr.evolution.io.Importer;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.NexusImporter.MissingBlockException;
import dr.evolution.io.NexusImporter.NexusBlock;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import org.jdom.JDOMException;

import java.io.*;
import java.util.*;
import java.util.List;


/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DataModelImporter {

    private final DateGuesser guesser;

    public DataModelImporter() {
        this.guesser = null;
    }

    public DataModelImporter(DateGuesser guesser) {
        this.guesser = guesser;
    }

    public HashMap importFromFile(File file) throws IOException, Importer.ImportException {
        HashMap dataModel = new HashMap();
        try {
            Reader reader = new FileReader(file);

            BufferedReader bufferedReader = new BufferedReader(reader);
            String line = bufferedReader.readLine();
            while (line != null && line.length() == 0) {
                line = bufferedReader.readLine();
            }

            if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
                // is a NEXUS file
                importNexusFile(file, guesser, dataModel);
            } else if ((line != null && line.trim().startsWith("" + FastaImporter.FASTA_FIRST_CHAR))) {
                // is a FASTA file
                importFastaFile(file, guesser, dataModel);
            } else if ((line != null && (line.toUpperCase().contains("<?XML") || line.toUpperCase().contains("<BEAST")))) {
                // assume it is a BEAST XML file and see if that works...
                importBEASTFile(file, guesser, dataModel);
//            } else {
//                // assume it is a tab-delimited traits file and see if that works...
//                importTraits(file);
            } else {
                throw new ImportException("Unrecognized format for imported file.");
            }

            bufferedReader.close();
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

        return dataModel;

    }

    // xml
    private void importBEASTFile(File file, DateGuesser guesser, Map dataModel) throws IOException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            BeastImporter importer = new BeastImporter(reader);

            List<TaxonList> taxonLists = new ArrayList<TaxonList>();
            List<Alignment> alignments = new ArrayList<Alignment>();

            importer.importBEAST(taxonLists, alignments);

            TaxonList taxa = taxonLists.get(0);

            int count = 1;
            for (Alignment alignment : alignments) {
                String name = file.getName();
                if (alignment.getId() != null && alignment.getId().length() > 0) {
                    name = alignment.getId();
                } else {
                    if (alignments.size() > 1) {
                        name += count;
                    }
                }
                setData(dataModel, guesser, name, taxa, taxonLists, alignment, null, null, null);

                count++;
            }

            reader.close();
        } catch (JDOMException e) {
            throw new ImportException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

    }

    // nexus
    private void importNexusFile(File file, DateGuesser guesser, Map dataModel) throws IOException, ImportException {
        TaxonList taxa = null;
        SimpleAlignment alignment = null;
        List<Tree> trees = new ArrayList<Tree>();
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

                        dataModel.put("taxa", createTaxonList(taxa));

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


                    } else if (block == NexusApplicationImporter.ASSUMPTIONS_BLOCK) {

                        importer.parseAssumptionsBlock(charSets);

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

            reader.close();

            // Allow the user to load taxa only (perhaps from a tree file) so that they can sample from a prior...
            if (alignment == null && taxa == null) {
                throw new MissingBlockException("TAXON, DATA or CHARACTERS block is missing");
            }

        } catch (IOException e) {
            throw new IOException(e.getMessage());
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
//        } catch (Exception e) {
//            throw new Exception(e.getMessage());
        }

        setData(dataModel, guesser, file.getName(), taxa, null, alignment, charSets, null, trees);
    }

    // FASTA

    private void importFastaFile(File file, DateGuesser guesser, Map dataModel) throws IOException, ImportException {
        try {
            FileReader reader = new FileReader(file);

            FastaImporter importer = new FastaImporter(reader, Nucleotides.INSTANCE);

            Alignment alignment = importer.importAlignment();

            reader.close();

            setData(dataModel, guesser, file.getName(), alignment, null, alignment, null, null, null);
        } catch (ImportException e) {
            throw new ImportException(e.getMessage());
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }
    }

    public Map importFromTreeFile(String fileName, Map dataModel) throws IOException, Importer.ImportException {
        Tree tree = null;
        try {
            Reader reader = new FileReader(fileName);

            BufferedReader bufferedReader = new BufferedReader(reader);
            String line = bufferedReader.readLine();
            while (line != null && line.length() == 0) {
                line = bufferedReader.readLine();
            }

            reader = new FileReader(fileName);

            if ((line != null && line.toUpperCase().contains("#NEXUS"))) {
                // is a NEXUS file
                NexusImporter importer = new NexusImporter(reader);
                tree = importer.importNextTree();

            } else {
                NewickImporter importer = new NewickImporter(reader);
                tree = importer.importNextTree();
            }

            bufferedReader.close();
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        }

        if (tree != null) {
            dataModel.put("tree", TreeUtils.newick(tree));
        }

        return dataModel;

    }


    private boolean isMissingValue(String value) {
        return (value.equals("?") || value.equals("NA") || value.length() == 0);
    }

//    public void importTraits(final File file, Map dataModel) throws Exception {
//        List<TraitData> importedTraits = new ArrayList<TraitData>();
//        Taxa taxa = options.taxonList;
//
//        DataTable<String[]> dataTable = DataTable.Text.parse(new FileReader(file));
//
//        String[] traitNames = dataTable.getColumnLabels();
//        String[] taxonNames = dataTable.getRowLabels();
//
//        for (int i = 0; i < dataTable.getColumnCount(); i++) {
//            boolean warningGiven = false;
//
//            String traitName = traitNames[i];
//
//            String[] values = dataTable.getColumn(i);
//            Class c = null;
//            if (!isMissingValue(values[0])) {
//                c = Utils.detectType(values[0]);
//            }
//            for (int j = 1; j < values.length; j++) {
//                if (!isMissingValue(values[j])) {
//                    if (c == null) {
//                        c = Utils.detectType(values[j]);
//                    } else {
//                        Class c1 = Utils.detectType(values[j]);
//                        if (c == Integer.class && c1 == Double.class) {
//                            // change the type to double
//                            c = Double.class;
//                        }
//
//                        if (c1 != c &&
//                                !(c == Double.class && c1 == Integer.class) &&
//                                !warningGiven ) {
//                            System.err.println("Not all values of same type for trait" + traitName);
//                            warningGiven = true;
//                        }
//                    }
//                }
//            }
//
//            TraitData.TraitType t = (c == Boolean.class || c == String.class || c == null) ? TraitData.TraitType.DISCRETE :
//                    (c == Integer.class) ? TraitData.TraitType.INTEGER : TraitData.TraitType.CONTINUOUS;
//            TraitData newTrait = new TraitData(options, traitName, file.getName(), t);
//
//            importedTraits.add(newTrait);
//
//            int j = 0;
//            for (final String taxonName : taxonNames) {
//
//                final int index = taxa.getTaxonIndex(taxonName);
//                Taxon taxon;
//                if (index >= 0) {
//                    taxon = taxa.getTaxon(index);
//                } else {
//                    taxon = new Taxon(taxonName);
//                    taxa.addTaxon(taxon);
//                }
//                if (!isMissingValue(values[j])) {
//                    taxon.setAttribute(traitName, Utils.constructFromString(c, values[j]));
//                } else {
//                    // AR - merge rather than replace existing trait values
//                    if (taxon.getAttribute(traitName) == null) {
//                        taxon.setAttribute(traitName, "?");
//                    }
//                }
//                j++;
//            }
//        }
//        setData(dataModel, file.getName(), taxa, null, null, null, importedTraits, null);
//    }

    // for Alignment
    private void setData(Map dataModel, DateGuesser guesser, String fileName, TaxonList taxonList, List<TaxonList> taxonLists, Alignment alignment,
                         List<NexusApplicationImporter.CharSet> charSets,
                         List<TraitData> traits, List<Tree> trees) throws ImportException, IllegalArgumentException {
        String fileNameStem = Utils.trimExtensions(fileName,
                new String[]{"NEX", "NEXUS", "FA", "FAS", "FASTA", "TRE", "TREE", "XML", "TXT"});

        checkTaxonList(taxonList, guesser);
        dataModel.put("taxa", createTaxonList(taxonList));
        dataModel.put("taxon_count", Integer.toString(taxonList.getTaxonCount()));
        dataModel.put("taxon_count_minus_one", Integer.toString(taxonList.getTaxonCount()-1));

        if (taxonLists != null) {
            List<Map> tss = new ArrayList<Map>();

            for (TaxonList tl : taxonLists) {
                Map ts = new HashMap();
                ts.put("id", tl.getId());
                ts.put("taxa", createTaxonList(taxonList));

                tss.add(ts);
            }
            dataModel.put("taxonSets", tss);
        }

        dataModel.put("alignment", createAlignment(alignment));
        dataModel.put("site_count", Integer.toString(alignment.getSiteCount()));

        dataModel.put("filename", fileName);
        dataModel.put("filename_stem", fileNameStem);
    }

    private void checkTaxonList(TaxonList taxonList, DateGuesser guesser) throws ImportException {

        // check the taxon names for invalid characters
        boolean foundAmp = false;
        for (Taxon taxon : taxonList) {
            String name = taxon.getId();
            if (name.indexOf('&') >= 0) {
                foundAmp = true;
            }
        }
        if (foundAmp) {
            throw new ImportException("One or more taxon names include an illegal character ('&').\n" +
                    "These characters will prevent BEAST from reading the resulting XML file.\n\n" +
                    "Please edit the taxon name(s) before reloading the data file.");
        }

        if (guesser != null) {
            guesser.guessDates(taxonList);
        } else {
            // make sure they all have dates...
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                if (taxonList.getTaxonAttribute(i, "date") == null) {
                    Date origin = new Date(0);

                    dr.evolution.util.Date date = dr.evolution.util.Date.createTimeSinceOrigin(0.0, Units.Type.YEARS, origin);
                    taxonList.getTaxon(i).setAttribute("date", date);
                }
            }
        }
    }

    private List<Map> createTaxonList(TaxonList taxa) {

        List<Map> tl = new ArrayList<Map>();
        for (Taxon taxon : taxa) {
            tl.add(createTaxon(taxon));
        }

        return tl;
    }

    private Map createTaxon(Taxon taxon) {
        Map t = new HashMap();
        t.put("id", taxon.getId());
        if (taxon.getDate() != null) {
            t.put("date", Double.toString(taxon.getDate().getTimeValue()));
        }
        return t;
    }

    private Map createAlignment(Alignment alignment) {

        Map a = new HashMap();
        a.put("id",  (alignment.getId() != null ? alignment.getId() : "alignment"));

        List<Map> ss = new ArrayList<Map>();
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            Sequence sequence = alignment.getSequence(i);
            ss.add(createSequence(sequence));
        }

        a.put("sequences", ss);
        return a;
    }

    private Map createSequence(Sequence sequence) {
        Map s = new HashMap();
        s.put("taxon", createTaxon(sequence.getTaxon()));
        s.put("data", sequence.getSequenceString());
        return s;
    }
}
