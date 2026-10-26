/*
 * BootstrapGenerator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.treelime.generator;

import dr.app.treelime.components.ComponentFactory;
import dr.app.treelime.options.AbstractPartitionData;
import dr.app.treelime.options.BootstrapReplicate;
import dr.app.treelime.options.PartitionData;
import dr.app.treelime.options.PartitionSubstitutionModel;
import dr.app.treelime.options.TreeLimeOptions;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SiteResampler;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Nonparametric bootstrap replicates of a TreeLime analysis written as files (plan 8 S3g): each replicate resamples
 * the columns of every partition with replacement, within each model partition (codon positions 1+2 pooled, or 1, 2
 * and 3 separately; padded constant sites included), and is written as an ordinary TreeLime XML whose partitions read
 * the resampled data. The starting tree still comes from the original data, so every replicate keeps the original
 * topology: a fixed-topology bootstrap. A counts file rebuilds any replicate.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BootstrapGenerator extends Generator {

    public static final String VERSION = "TreeLime bootstrap 1";
    public static final String COUNTS_FILE = "bootstrap.counts.tsv";
    public static final String SCRIPT_FILE = "run_replicates.sh";

    private final ComponentFactory[] components;

    public BootstrapGenerator(TreeLimeOptions options, ComponentFactory[] components) {
        super(options, components);
        this.components = components;
    }

    /**
     * One partition's draw: how often each of its sites was drawn and, when it is padded, the redrawn constant-site
     * counts.
     */
    public static final class PartitionDraw {
        public final PartitionData partition;
        public final int[] siteCounts;
        public final double[] constantCounts;

        PartitionDraw(PartitionData partition, int[] siteCounts, double[] constantCounts) {
            this.partition = partition;
            this.siteCounts = siteCounts;
            this.constantCounts = constantCounts;
        }
    }

    /**
     * @return the partitions the replicates resample, after checking that they can be
     */
    public List<PartitionData> getPartitions() throws GeneratorException {
        if (options.samplePriorOnly) {
            throw new GeneratorException("Bootstrap replicates need the data: sampling from the prior only " +
                    "is selected.");
        }
        List<PartitionData> partitions = new ArrayList<PartitionData>();
        for (AbstractPartitionData data : options.dataPartitions) {
            if (!(data instanceof PartitionData) || ((PartitionData) data).getAlignment() == null) {
                throw new GeneratorException("Bootstrap replicates resample sequence alignments only; " +
                        "partition " + data.getName() + " is not one.");
            }
            int type = data.getDataType().getType();
            if (type != DataType.NUCLEOTIDES && type != DataType.AMINO_ACIDS) {
                throw new GeneratorException("Bootstrap replicates support nucleotide and amino-acid " +
                        "partitions; partition " + data.getName() + " is " + data.getDataType().getDescription() + ".");
            }
            partitions.add((PartitionData) data);
        }
        if (partitions.isEmpty()) {
            throw new GeneratorException("There are no alignments to resample.");
        }
        return partitions;
    }

    /**
     * @return the 0-based alignment columns of a partition, in order
     */
    public static int[] sites(PartitionData partition) {
        Alignment alignment = partition.getAlignment();
        int from = Math.max(1, partition.getFromSite()) - 1;
        int to = partition.getToSite() < 1 ? alignment.getSiteCount() : partition.getToSite();
        int every = Math.max(1, partition.getEvery());
        int[] sites = new int[(to - from + every - 1) / every];
        for (int i = 0; i < sites.length; i++) {
            sites[i] = from + i * every;
        }
        return sites;
    }

    /**
     * @return the model partitions of a partition, as indices into its sites: one group, or the codon positions
     * (1+2 pooled under the 112 pattern)
     */
    public static List<int[]> groups(PartitionData partition) {
        int n = sites(partition).length;
        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
        String pattern = model.getCodonHeteroPattern();
        List<int[]> groups = new ArrayList<int[]>();
        if (partition.getDataType().getType() == DataType.NUCLEOTIDES && pattern != null &&
                model.getCodonPartitionCount() > 1) {
            int[][] positions = pattern.equals("112") ? new int[][]{{0, 1}, {2}} : new int[][]{{0}, {1}, {2}};
            for (int[] codonPositions : positions) {
                List<Integer> group = new ArrayList<Integer>();
                for (int i = 0; i < n; i++) {
                    for (int position : codonPositions) {
                        if (i % 3 == position) {
                            group.add(i);
                        }
                    }
                }
                int[] indices = new int[group.size()];
                for (int i = 0; i < indices.length; i++) {
                    indices[i] = group.get(i);
                }
                groups.add(indices);
            }
        } else {
            int[] all = new int[n];
            for (int i = 0; i < n; i++) {
                all[i] = i;
            }
            groups.add(all);
        }
        return groups;
    }

    /**
     * @return replicate b's draw for every partition, from the master seed alone
     */
    public List<PartitionDraw> draw(int replicate, long masterSeed) throws GeneratorException {
        Random random = new Random(SiteResampler.replicateSeed(masterSeed, replicate));
        List<PartitionDraw> draws = new ArrayList<PartitionDraw>();
        for (PartitionData partition : getPartitions()) {
            int n = sites(partition).length;
            int[] siteCounts = new int[n];
            double[] padding = partition.getConstantSiteCounts();
            double[] constantCounts = null;
            if (padding != null) {
                // the padded constant sites are columns too: draw from the sites and the constants together
                int constant = 0;
                for (double count : padding) {
                    constant += (int) count;
                }
                int[] counts = SiteResampler.resample(n + constant, random);
                System.arraycopy(counts, 0, siteCounts, 0, n);
                constantCounts = new double[padding.length];
                int column = n;
                for (int state = 0; state < padding.length; state++) {
                    for (int k = 0; k < (int) padding[state]; k++) {
                        constantCounts[state] += counts[column++];
                    }
                }
            } else {
                for (int[] group : groups(partition)) {
                    int[] counts = SiteResampler.resample(group.length, random);
                    for (int k = 0; k < group.length; k++) {
                        siteCounts[group[k]] = counts[k];
                    }
                }
            }
            draws.add(new PartitionDraw(partition, siteCounts, constantCounts));
        }
        return draws;
    }

    /**
     * @return the draw in which every site is taken once: the original data
     */
    public List<PartitionDraw> identity() throws GeneratorException {
        List<PartitionDraw> draws = new ArrayList<PartitionDraw>();
        for (PartitionData partition : getPartitions()) {
            int[] siteCounts = new int[sites(partition).length];
            Arrays.fill(siteCounts, 1);
            draws.add(new PartitionDraw(partition, siteCounts, partition.getConstantSiteCounts()));
        }
        return draws;
    }

    /**
     * The resampled alignment of one partition. Within each model partition the drawn columns fill that group's
     * positions in order, each column as often as it was drawn, so the counts alone rebuild the alignment and codon
     * positions stay where the patterns read them. A padded partition is a single group of variable length.
     */
    public static Alignment replicateAlignment(PartitionDraw draw) {
        PartitionData partition = draw.partition;
        Alignment original = partition.getAlignment();
        int[] sites = sites(partition);
        List<Integer> columns = new ArrayList<Integer>();
        if (draw.constantCounts != null) {
            for (int i = 0; i < sites.length; i++) {
                for (int k = 0; k < draw.siteCounts[i]; k++) {
                    columns.add(sites[i]);
                }
            }
        } else {
            Integer[] slots = new Integer[sites.length];
            for (int[] group : groups(partition)) {
                int slot = 0;
                for (int k = 0; k < group.length; k++) {
                    for (int c = 0; c < draw.siteCounts[group[k]]; c++) {
                        slots[group[slot++]] = sites[group[k]];
                    }
                }
            }
            columns.addAll(Arrays.asList(slots));
        }

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(original.getDataType());
        for (int t = 0; t < original.getSequenceCount(); t++) {
            String sequence = original.getAlignedSequenceString(t);
            StringBuilder resampled = new StringBuilder(columns.size());
            for (int column : columns) {
                resampled.append(sequence.charAt(column));
            }
            Sequence replicate = new Sequence(original.getTaxon(t), resampled.toString());
            replicate.setDataType(original.getDataType());
            alignment.addSequence(replicate);
        }
        alignment.setId(BootstrapReplicate.alignmentId(partition));
        return alignment;
    }

    /**
     * Writes replicate b's XML into the directory; the same master seed always gives the same file.
     */
    public File writeReplicate(File directory, int replicate, long masterSeed) throws Exception {
        prepare();
        return writeReplicate(directory, replicate, draw(replicate, masterSeed));
    }

    // as before writing any TreeLime XML: bring the parameter and operator lists up to date, then check the options
    private void prepare() throws GeneratorException {
        options.selectParameters();
        options.selectOperators();
        new BeastGenerator(options, components).checkOptions();
    }

    public File writeReplicate(File directory, int replicate, List<PartitionDraw> draws) throws Exception {
        BootstrapReplicate data = new BootstrapReplicate();
        for (PartitionDraw draw : draws) {
            data.put(draw.partition, replicateAlignment(draw), draw.constantCounts);
        }
        String stem = options.fileNameStem;
        String reportFileName = options.maximizationOptions.reportFileName;
        File file = new File(directory, "rep" + replicate + ".xml");
        try {
            options.bootstrapReplicate = data;
            // the outputs of each replicate carry its number
            options.fileNameStem = stem + ".rep" + replicate;
            options.maximizationOptions.reportFileName = null;
            new BeastGenerator(options, components).generateXML(file);
        } finally {
            options.bootstrapReplicate = null;
            options.fileNameStem = stem;
            options.maximizationOptions.reportFileName = reportFileName;
        }
        return file;
    }

    /**
     * Writes replicates 1 to B, the counts file and the run script.
     */
    public void generate(File directory, int replicates, long masterSeed) throws Exception {
        if (replicates < 1) {
            throw new GeneratorException("Choose at least one bootstrap replicate.");
        }
        List<PartitionData> partitions = getPartitions();
        prepare();
        if (!directory.isDirectory() && !directory.mkdirs()) {
            throw new IOException("Cannot create the directory " + directory);
        }

        PrintWriter counts = new PrintWriter(new FileWriter(new File(directory, COUNTS_FILE)));
        try {
            counts.println("# " + VERSION + ": how often each column of each partition was drawn");
            counts.println("# master_seed=" + masterSeed + "\treplicates=" + replicates + "\tstem=" + options.fileNameStem);
            counts.println("# fixed-topology bootstrap: every replicate keeps the original starting tree");
            for (PartitionData partition : partitions) {
                int[] sites = sites(partition);
                counts.println("# partition " + partition.getName() + ": alignment columns " + (sites[0] + 1) + "-" +
                        (sites[sites.length - 1] + 1) + " every " + Math.max(1, partition.getEvery()) +
                        (groups(partition).size() > 1 ? ", codon positions " +
                                partition.getPartitionSubstitutionModel().getCodonHeteroPattern() : "") +
                        (partition.getConstantSiteCounts() != null ? ", padded with constant sites" : ""));
            }
            counts.println("replicate\tpartition\tcolumn_counts\tconstant_counts");
            for (int b = 1; b <= replicates; b++) {
                List<PartitionDraw> draws = draw(b, masterSeed);
                writeReplicate(directory, b, draws);
                for (PartitionDraw draw : draws) {
                    StringBuilder row = new StringBuilder().append(b).append('\t').append(draw.partition.getName())
                            .append('\t');
                    for (int i = 0; i < draw.siteCounts.length; i++) {
                        row.append(i > 0 ? "," : "").append(draw.siteCounts[i]);
                    }
                    row.append('\t');
                    if (draw.constantCounts != null) {
                        for (int i = 0; i < draw.constantCounts.length; i++) {
                            row.append(i > 0 ? "," : "").append((long) draw.constantCounts[i]);
                        }
                    }
                    counts.println(row);
                }
            }
        } finally {
            counts.close();
        }
        writeScript(directory, replicates);
    }

    private void writeScript(File directory, int replicates) throws IOException {
        File script = new File(directory, SCRIPT_FILE);
        PrintWriter out = new PrintWriter(new FileWriter(script));
        try {
            out.println("#!/bin/sh");
            out.println("# Runs the " + replicates + " bootstrap replicates of " + options.fileNameStem +
                    " with BEAST, one after another.");
            out.println("# Fixed-topology bootstrap: every replicate keeps the original starting tree.");
            out.println("# Set BEAST to the BEAST command if it is not 'beast', e.g.");
            out.println("#   BEAST=\"java -Djava.library.path=/usr/local/lib -jar /path/to/beast.jar\" sh " + SCRIPT_FILE);
            out.println("cd \"$(dirname \"$0\")\" || exit 1");
            out.println("BEAST=\"${BEAST:-beast}\"");
            out.println("for b in $(seq 1 " + replicates + "); do");
            out.println("    $BEAST -overwrite -seed \"$b\" \"rep$b.xml\" > \"rep$b.log\" 2>&1 || " +
                    "echo \"replicate $b failed: see rep$b.log\"");
            out.println("done");
        } finally {
            out.close();
        }
        script.setExecutable(true);
    }
}
