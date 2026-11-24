/*
 * Copyright (c) 2026 the BEAST Development Team.
 * Licensed under the GNU Lesser General Public License, version 2 or later.
 */
package dr.app.tools;

import beagle.Beagle;
import beagle.BeagleFlag;
import beagle.InstanceDetails;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.HierarchicalSubstitutionModelParser;
import dr.inference.model.Parameter;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** Ordinary BEAST likelihood evaluations with standard, separate per-branch buffers. */
public final class HierarchicalPracticalBenchmark {
    private enum Method { H, HA, NE, AE }

    private static final class Options {
        Path input, output;
        Method method;
        String endpoint = "model-update";
        int warmupEvals = 20, batches = 5, blockSize = 0;
        long maxEvals = 100000000L;
        double warmupSeconds = 1.0, batchSeconds = 0.5;
        boolean validateOnly;

        static Options parse(String[] args) {
            Options o = new Options();
            for (int i = 0; i < args.length; i++) {
                String key = args[i];
                if ("--validate-only".equals(key)) { o.validateOnly = true; continue; }
                if (i + 1 == args.length) throw new IllegalArgumentException("Missing value for " + key);
                String value = args[++i];
                switch (key) {
                    case "--input": o.input = Paths.get(value).toAbsolutePath(); break;
                    case "--output": o.output = Paths.get(value).toAbsolutePath(); break;
                    case "--method": o.method = Method.valueOf(value); break;
                    case "--endpoint": o.endpoint = value; break;
                    case "--warmup-evals": o.warmupEvals = Integer.parseInt(value); break;
                    case "--warmup-seconds": o.warmupSeconds = Double.parseDouble(value); break;
                    case "--batch-seconds": o.batchSeconds = Double.parseDouble(value); break;
                    case "--batches": o.batches = Integer.parseInt(value); break;
                    case "--pattern-block-size": o.blockSize = Integer.parseInt(value); break;
                    case "--max-evals": o.maxEvals = Long.parseLong(value); break;
                    default: throw new IllegalArgumentException("Unknown option " + key);
                }
            }
            if (o.input == null || o.output == null || o.method == null)
                throw new IllegalArgumentException("Required: --input FILE --method H|HA|NE|AE --output FILE");
            if (!o.endpoint.equals("model-update"))
                throw new IllegalArgumentException("--endpoint must be model-update for the practical benchmark");
            if (!Double.isFinite(o.warmupSeconds) || !Double.isFinite(o.batchSeconds)
                    || o.warmupEvals < 0 || o.warmupSeconds < 0 || o.batchSeconds <= 0 || o.batches < 1
                    || o.blockSize < 0 || o.maxEvals < 1) throw new IllegalArgumentException("Invalid timing/block limit");
            return o;
        }
    }

    private static final class Input {
        int s, n, l, k, root;
        int[] left, right, operations, matrices;
        double[] rawRates, pi, branches, matrixLengths, weights;
        int[][] tips;
        boolean sharedTransitions;
        int distinctBranchLengths;

        static Input read(Path path) throws IOException {
            Input x = new Input();
            try (DataInputStream in = new DataInputStream(new BufferedInputStream(Files.newInputStream(path)))) {
                byte[] magic = new byte[8]; in.readFully(magic);
                if (!Arrays.equals(magic, "JHRATE01".getBytes(StandardCharsets.US_ASCII)))
                    throw new IOException("Not a JHRATE01 input");
                x.s = in.readInt(); x.n = in.readInt(); x.l = in.readInt(); x.k = in.readInt(); x.root = in.readInt();
                if (x.s < 2 || x.n < 2 || x.l < 1 || x.k < 1 || x.k > x.l || x.root != 2 * x.s - 2)
                    throw new IOException("Invalid input dimensions/root");
                x.pi = new double[x.s];
                for (int i = 0; i < x.s; i++) x.pi[i] = in.readDouble();
                x.left = new int[2 * x.s - 1]; x.right = new int[x.left.length];
                x.rawRates = new double[x.left.length];
                for (int i = x.s; i < x.left.length; i++) {
                    int parent = in.readInt(), left = in.readInt(), right = in.readInt();
                    double rate = in.readDouble();
                    if (parent != i || left < 0 || right < 0 || left >= i || right >= i || left == right
                            || !Double.isFinite(rate) || rate < 0) throw new IOException("Invalid character record " + i);
                    x.left[i] = left; x.right[i] = right; x.rawRates[i] = rate;
                }
                x.operations = new int[7 * (x.n - 1)];
                for (int i = 0; i < x.n - 1; i++) {
                    int parent = in.readInt(), left = in.readInt(), right = in.readInt();
                    if (parent != x.n + i || left < 0 || right < 0 || left >= parent || right >= parent || left == right)
                        throw new IOException("Invalid phylogeny record " + i);
                    int p = 7 * i;
                    x.operations[p] = parent; x.operations[p + 1] = x.operations[p + 2] = Beagle.NONE;
                    x.operations[p + 3] = left; x.operations[p + 4] = left;
                    x.operations[p + 5] = right; x.operations[p + 6] = right;
                }
                x.branches = new double[2 * x.n - 2]; x.matrices = new int[x.branches.length];
                for (int i = 0; i < x.branches.length; i++) {
                    x.branches[i] = in.readDouble(); x.matrices[i] = i;
                    if (!(x.branches[i] >= 0) || !Double.isFinite(x.branches[i])) throw new IOException("Invalid branch length");
                }
                Set<Long> distinctLengths = new HashSet<>();
                for (double branch : x.branches) distinctLengths.add(Double.doubleToLongBits(branch));
                x.distinctBranchLengths = distinctLengths.size();
                // Ordinary analysis behavior: every edge owns its matrix even when lengths happen to match.
                x.sharedTransitions = false;
                x.matrixLengths = x.branches;
                x.weights = new double[x.k]; double sum = 0;
                for (int i = 0; i < x.k; i++) {
                    x.weights[i] = in.readDouble(); sum += x.weights[i];
                    if (!(x.weights[i] > 0) || x.weights[i] != Math.rint(x.weights[i])) throw new IOException("Invalid pattern weight");
                }
                if (sum != x.l) throw new IOException("Pattern weights do not sum to sites");
                x.tips = new int[x.n][x.k];
                for (int i = 0; i < x.n; i++) for (int j = 0; j < x.k; j++) {
                    int state = in.readInt();
                    if (state < 0 || state >= x.s) throw new IOException("Invalid observed state");
                    x.tips[i][j] = state;
                }
                if (in.read() != -1) throw new IOException("Trailing input bytes");
            }
            return x;
        }

        String newick(int node) {
            return node < s ? String.format(Locale.ROOT, "state%03d:1", node)
                    : "(" + newick(left[node]) + "," + newick(right[node]) + "):1";
        }

        Map<BitSet, Double> ratesByClade() {
            Map<BitSet, Double> result = new HashMap<>();
            BitSet[] clades = new BitSet[2 * s - 1];
            for (int i = 0; i < s; i++) { clades[i] = new BitSet(s); clades[i].set(i); }
            for (int i = s; i < clades.length; i++) {
                clades[i] = (BitSet) clades[left[i]].clone(); clades[i].or(clades[right[i]]);
                result.put(clades[i], rawRates[i]);
            }
            return result;
        }
    }


    /** Import a verified table of unique weighted patterns without changing their values. */
    private static final class WeightedPatterns extends Patterns {
        WeightedPatterns(Input input, GeneralDataType type, Taxa taxa) {
            super(type, taxa);
            patternCount = input.k;
            patterns = new int[input.k][input.n];
            weights = input.weights.clone();
            for (int k = 0; k < input.k; k++)
                for (int tip = 0; tip < input.n; tip++) patterns[k][tip] = input.tips[tip][k];
            setId("patterns");
        }
    }

    private static final class PatternKey {
        final int[] states;
        final int hash;
        PatternKey(int[] states) { this.states = states; hash = Arrays.hashCode(states); }
        @Override public int hashCode() { return hash; }
        @Override public boolean equals(Object other) {
            return other instanceof PatternKey && Arrays.equals(states, ((PatternKey) other).states);
        }
    }

    private final Options options;
    private final Input input;
    private final Parameter rates;
    private final int rootRateIndex;
    private final HierarchicalSubstitutionModel model;
    private final TreeModel phylogeny;
    private final BeagleDataLikelihoodDelegate delegate;
    private final TreeDataLikelihood likelihood;
    private final int[] inputPatternOf;
    private int point;

    private HierarchicalPracticalBenchmark(Options options, Input input) throws Exception {
        this.options = options; this.input = input;
        if (options.blockSize > 0 && options.blockSize < input.k)
            throw new IllegalArgumentException("Practical BEAST benchmark uses the complete resident pattern list; streaming is diagnostic only");
        String[] labels = new String[input.s];
        for (int i = 0; i < input.s; i++) labels[i] = String.format(Locale.ROOT, "state%03d", i);
        GeneralDataType dataType = new GeneralDataType(labels);
        TreeModel characterTree = new DefaultTreeModel("characterTree",
                new NewickImporter(input.newick(input.root) + ";").importTree(null));
        rates = new Parameter.Default(input.s - 1, 1.0);
        fillRates(characterTree, characterTree.getRoot(), input.ratesByClade());
        rootRateIndex = characterTree.getRoot().getNumber() - input.s;
        FrequencyModel frequency = new FrequencyModel(dataType, new Parameter.Default(input.pi.clone()));
        model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, characterTree, frequency, true));
        model.setEigenSystemType(options.method == Method.NE
                ? HierarchicalSubstitutionModel.EigenSystemType.NUMERICAL : HierarchicalSubstitutionModel.EigenSystemType.ANALYTICAL);
        model.setAggregation(options.method == Method.HA ? HierarchicalSubstitutionModel.Aggregation.ON
                : HierarchicalSubstitutionModel.Aggregation.OFF);
        boolean hierarchical = options.method == Method.H || options.method == Method.HA;
        long required = BeagleFlag.PROCESSOR_CPU.getMask() | BeagleFlag.PRECISION_DOUBLE.getMask()
                | BeagleFlag.THREADING_NONE.getMask();
        if (hierarchical) required |= BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()
                | (options.method == Method.HA ? BeagleFlag.HIERARCHICAL_AGGREGATE_ON.getMask()
                : BeagleFlag.HIERARCHICAL_AGGREGATE_OFF.getMask());
        else required |= BeagleFlag.EIGEN_REAL.getMask() | BeagleFlag.VECTOR_SSE.getMask();
        System.setProperty("beagle.preferred.flags", Long.toString(required));
        System.setProperty("beagle.required.flags", Long.toString(required));
        System.setProperty("beagle.thread.count", "0");
        System.setProperty("beagle.force.vectorization", "true");
        System.setProperty("beagle.resource.auto", "false");
        phylogeny = new DefaultTreeModel("phylogeny", new NewickImporter(phylogenyNewick(2 * input.n - 2) + ";").importTree(null));
        Taxa taxa = new Taxa();
        for (int i = 0; i < input.n; i++) taxa.addTaxon(new Taxon(String.format(Locale.ROOT, "taxon%03d", i)));
        PatternList original = new WeightedPatterns(input, dataType, taxa);
        HomogeneousBranchModel branchModel = new HomogeneousBranchModel(model);
        PatternList processed = hierarchical
                ? HierarchicalSubstitutionModelParser.aggregatedPatterns(original, branchModel) : original;
        if (processed.getPatternCount() != input.k)
            throw new IllegalStateException("Normal pattern processing changed the unique pattern count");
        inputPatternOf = new int[input.k];
        Map<PatternKey, Integer> index = new HashMap<>();
        for (int k = 0; k < input.k; k++) {
            if (index.put(new PatternKey(original.getPattern(k)), k) != null)
                throw new IllegalArgumentException("Input patterns must be unique");
        }
        for (int k = 0; k < processed.getPatternCount(); k++) {
            Integer originalIndex = index.get(new PatternKey(processed.getPattern(k)));
            if (originalIndex == null || processed.getPatternWeight(k) != input.weights[originalIndex])
                throw new IllegalStateException("Normal pattern processing changed a pattern or weight");
            inputPatternOf[k] = originalIndex;
        }
        delegate = new BeagleDataLikelihoodDelegate(phylogeny, processed, branchModel,
                new GammaSiteRateModel("siteRate"), false, false, PartialsRescalingScheme.DEFAULT, true,
                new PreOrderSettings(false, false, false, false, false, false, false, hierarchical));
        likelihood = new TreeDataLikelihood(delegate, phylogeny, new DefaultBranchRateModel());
        InstanceDetails details = delegate.getBeagleInstance().getDetails();
        if (details == null || (details.getFlags() & required) != required)
            throw new IllegalStateException("Required native implementation unavailable");
        if (!hierarchical && (details.getFlags() & BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()) != 0)
            throw new IllegalStateException("Dense method selected hierarchical representation");
        Set<Integer> matrices = branchMatrixIndices();
        if (matrices.size() != 2 * input.n - 2)
            throw new IllegalStateException("Ordinary BEAST delegate did not assign a distinct matrix to every branch");
    }

    private Set<Integer> branchMatrixIndices() {
        Set<Integer> matrices = new TreeSet<>();
        EvolutionaryProcessDelegate process = delegate.getEvolutionaryProcessDelegate();
        for (int i = 0; i < phylogeny.getNodeCount(); i++)
            if (!phylogeny.isRoot(phylogeny.getNode(i))) matrices.add(process.getMatrixIndex(i));
        return matrices;
    }

    private String phylogenyNewick(int node) {
        String value;
        if (node < input.n) value = String.format(Locale.ROOT, "taxon%03d", node);
        else {
            int op = 7 * (node - input.n);
            value = "(" + phylogenyNewick(input.operations[op + 3]) + ","
                    + phylogenyNewick(input.operations[op + 5]) + ")";
        }
        return node == 2 * input.n - 2 ? value : value + ":" + Double.toString(input.branches[node]);
    }

    private BitSet fillRates(TreeModel tree, NodeRef node, Map<BitSet, Double> byClade) {
        BitSet clade = new BitSet(input.s);
        if (tree.isExternal(node)) clade.set(Integer.parseInt(tree.getTaxonId(node.getNumber()).substring(5)));
        else {
            for (int c = 0; c < tree.getChildCount(node); c++) clade.or(fillRates(tree, tree.getChild(node, c), byClade));
            Double rate = byClade.get(clade);
            if (rate == null) throw new IllegalArgumentException("Character rate mapping mismatch");
            rates.setParameterValue(node.getNumber() - input.s, rate);
        }
        return clade;
    }

    private double evaluate(boolean update) {
        if (update) {
            point = 1 - point;
            rates.setParameterValue(rootRateIndex, input.rawRates[input.root] * (point == 0 ? 1.0 : 1.01));
        }
        double result = likelihood.getLogLikelihood();
        if (!Double.isFinite(result)) throw new ArithmeticException("Nonfinite practical BEAST likelihood");
        return result;
    }

    private double[] siteLogs() {
        double[] actual = delegate.getSiteLogLikelihoods(), restored = new double[input.k];
        if (actual.length != input.k) throw new IllegalStateException("Pattern count changed");
        for (int k = 0; k < actual.length; k++) restored[inputPatternOf[k]] = actual[k];
        return restored;
    }

    private static void saveDoubles(Path path, double[] values) throws IOException {
        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(path)))) {
            for (double value : values) {
                if (!Double.isFinite(value)) throw new ArithmeticException("Nonfinite pattern likelihood");
                out.writeDouble(value);
            }
        }
    }

    private long matrixUpdates() {
        Matcher match = Pattern.compile("matrix updates = (\\d+)").matcher(delegate.getReport());
        if (!match.find()) throw new IllegalStateException("Normal delegate counters unavailable");
        return Long.parseLong(match.group(1));
    }

    private Map<String, Object> run(long mainStart) throws IOException {
        Map<String, Object> result = map();
        double base = evaluate(false);
        double preparation = (System.nanoTime() - mainStart) / 1e9;
        Map<String, Object> marker = map(); marker.put("log_likelihood", base); marker.put("java_main_to_first_seconds", preparation);
        System.out.println("FIRST_LIKELIHOOD " + json(marker)); System.out.flush();
        Path baseFile = options.output.resolveSibling(options.output.getFileName() + ".base.f64be");
        Path changedFile = options.output.resolveSibling(options.output.getFileName() + ".changed.f64be");
        saveDoubles(baseFile, siteLogs());
        double changed = evaluate(true); saveDoubles(changedFile, siteLogs());
        evaluate(true);
        Map<String, Object> validation = map();
        validation.put("base_log_likelihood", base); validation.put("changed_log_likelihood", changed);
        validation.put("base_site_logs", baseFile.toString()); validation.put("changed_site_logs", changedFile.toString());
        validation.put("site_log_dtype", ">f8"); result.put("validation", validation);
        result.put("method", options.method.name()); result.put("endpoint", options.endpoint);
        result.put("benchmark_backend", "ordinary BEAST TreeDataLikelihood and BeagleDataLikelihoodDelegate");
        result.put("execution_mode", "resident"); result.put("pattern_blocks", 1); result.put("padding_patterns", 0);
        result.put("states", input.s); result.put("tips", input.n); result.put("sites", input.l); result.put("patterns", input.k);
        result.put("preparation_seconds", preparation);
        result.put("preparation_clock_origin", "Java main entry; FIRST_LIKELIHOOD marker permits external process launch bound");
        result.put("first_log_likelihood", base); result.put("native_patterns", input.k);
        result.put("shared_transition_matrices", false); result.put("branch_count", 2 * input.n - 2);
        result.put("matrix_buffer_count", delegate.getEvolutionaryProcessDelegate().getMatrixBufferCount());
        result.put("distinct_branch_matrix_indices", branchMatrixIndices());
        result.put("partial_buffer_count", delegate.getPartialBufferCount());
        result.put("scaling", delegate.getRescalingScheme().toString());
        result.put("delay_rescaling_until_underflow", delegate.getDelayRescalingUntilUnderflow());
        result.put("underflow_count_after_validation", delegate.getUnderflowCount());
        result.put("category_count", 1); result.put("threads", 1);
        result.put("root_rate_multiplier_base", 1.0); result.put("root_rate_multiplier_changed", 1.01);
        result.put("input", options.input.toString()); result.put("input_sha256", sha256(options.input));
        InstanceDetails details = delegate.getBeagleInstance().getDetails();
        boolean hierarchical = options.method == Method.H || options.method == Method.HA;
        result.put("implementation", details.getImplementationName() == null
                ? (hierarchical ? "CPU-Hierarchical-Double" : input.s == 4 ? "CPU-4State-SSE-Double" : "CPU-SSE-Double")
                : details.getImplementationName());
        result.put("implementation_name_source", details.getImplementationName() == null ? "strict native flags" : "native JNI");
        result.put("flags", details.getFlags()); result.put("flag_names", BeagleFlag.toString(details.getFlags()));
        result.put("resource", details.getResourceNumber());
        result.put("aggregation_requested", options.method == Method.HA ? "on" : options.method == Method.H ? "off" : "not applicable");
        result.put("pattern_order", "normal BEAST pattern processing; per-pattern outputs restored to input order");
        result.put("transition_reuse_note", "Ordinary BEAST per-branch matrices and standard buffer management; no equal-length sharing");
        List<Object> batches = new ArrayList<>(); result.put("batches", batches);
        if (!options.validateOnly) {
            long start = System.nanoTime(), count = 0;
            do {
                evaluate(true); count++;
                if (count >= options.maxEvals && (count < options.warmupEvals || (System.nanoTime() - start) / 1e9 < options.warmupSeconds))
                    throw new IllegalStateException("Warmup max-evals reached before minimum");
            } while (count < options.warmupEvals || (System.nanoTime() - start) / 1e9 < options.warmupSeconds);
            Map<String, Object> warmup = map(); warmup.put("evaluations", count); warmup.put("seconds", (System.nanoTime() - start) / 1e9);
            result.put("warmup", warmup);
            double[] times = new double[options.batches];
            for (int b = 0; b < options.batches; b++) {
                long matrixBefore = matrixUpdates(), partialsBefore = delegate.getTotalCalculationCount();
                start = System.nanoTime(); count = 0; double last;
                do {
                    last = evaluate(true); count++;
                    if (count >= options.maxEvals && (System.nanoTime() - start) / 1e9 < options.batchSeconds)
                        throw new IllegalStateException("Batch max-evals reached before minimum");
                } while ((System.nanoTime() - start) / 1e9 < options.batchSeconds || (count & 1) != 0);
                double seconds = (System.nanoTime() - start) / 1e9;
                long matrices = matrixUpdates() - matrixBefore;
                long partials = delegate.getTotalCalculationCount() - partialsBefore;
                if (matrices < count * (2 * input.n - 2) || partials < count * (input.n - 1))
                    throw new IllegalStateException("A practical evaluation was cached or skipped branch updates");
                Map<String, Object> batch = map(); batch.put("batch", b + 1); batch.put("iterations", count);
                batch.put("seconds", seconds); batch.put("ms_per_eval", seconds * 1000.0 / count);
                batch.put("last_log_likelihood", last); batch.put("updated_transition_matrices", matrices);
                batch.put("phylogeny_operations", partials); batches.add(batch);
                times[b] = seconds * 1000.0 / count;
            }
            Arrays.sort(times); result.put("median_ms_per_eval", times.length % 2 == 1 ? times[times.length / 2]
                    : 0.5 * (times[times.length / 2 - 1] + times[times.length / 2]));
        }
        result.put("normal_delegate_report", delegate.getReport());
        result.put("underflow_count", delegate.getUnderflowCount());
        result.put("status", "ok"); result.put("format_version", 2);
        return result;
    }

    private static Map<String, Object> map() { return new LinkedHashMap<>(); }

    private static String sha256(Path path) throws IOException {
        try {
            MessageDigest hash = MessageDigest.getInstance("SHA-256");
            byte[] buffer = new byte[1024 * 1024];
            try (InputStream in = Files.newInputStream(path)) { int n; while ((n = in.read(buffer)) >= 0) hash.update(buffer, 0, n); }
            StringBuilder result = new StringBuilder();
            for (byte b : hash.digest()) result.append(String.format(Locale.ROOT, "%02x", b & 255));
            return result.toString();
        } catch (java.security.NoSuchAlgorithmException e) { throw new AssertionError(e); }
    }

    private static String json(Object value) {
        if (value == null) return "null";
        if (value instanceof Number || value instanceof Boolean) return value.toString();
        if (value instanceof Map) {
            StringBuilder out = new StringBuilder("{"); boolean first = true;
            for (Object entryObject : ((Map<?, ?>) value).entrySet()) {
                Map.Entry<?, ?> entry = (Map.Entry<?, ?>) entryObject;
                if (!first) out.append(','); first = false;
                out.append(json(entry.getKey().toString())).append(':').append(json(entry.getValue()));
            }
            return out.append('}').toString();
        }
        if (value instanceof Iterable) {
            StringBuilder out = new StringBuilder("["); boolean first = true;
            for (Object item : (Iterable<?>) value) { if (!first) out.append(','); first = false; out.append(json(item)); }
            return out.append(']').toString();
        }
        String text = value.toString(); StringBuilder out = new StringBuilder("\"");
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (c == '\\' || c == '"') out.append('\\').append(c);
            else if (c == '\n') out.append("\\n"); else if (c == '\r') out.append("\\r");
            else if (c == '\t') out.append("\\t"); else if (c < 32) out.append(String.format(Locale.ROOT, "\\u%04x", (int) c));
            else out.append(c);
        }
        return out.append('"').toString();
    }

    private static void saveJson(Path path, Map<String, Object> value) throws IOException {
        Path temporary = path.resolveSibling(path.getFileName() + ".tmp");
        Files.write(temporary, (json(value) + "\n").getBytes(StandardCharsets.UTF_8));
        Files.move(temporary, path, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE);
    }

    public static void main(String[] args) throws Throwable {
        long mainStart = System.nanoTime(); Options options = null; HierarchicalPracticalBenchmark benchmark = null;
        try {
            options = Options.parse(args); Files.createDirectories(options.output.getParent());
            if (Boolean.getBoolean("java.only")) throw new IllegalArgumentException("Native BEAGLE is required");
            if (System.getenv("BEAGLE_HRATE_AGGREGATE") != null)
                throw new IllegalArgumentException("Unset BEAGLE_HRATE_AGGREGATE; ordinary model configuration selects aggregation");
            benchmark = new HierarchicalPracticalBenchmark(options, Input.read(options.input));
            Map<String, Object> result = benchmark.run(mainStart); saveJson(options.output, result);
            System.out.println(json(result));
        } catch (Throwable error) {
            if (options != null && options.output != null) {
                Map<String, Object> failure = map(); failure.put("status", "failed");
                failure.put("method", options.method == null ? null : options.method.name());
                failure.put("endpoint", options.endpoint); failure.put("error", error.toString());
                saveJson(options.output, failure);
            }
            throw error;
        } finally {
            if (benchmark != null) BeagleDataLikelihoodDelegate.releaseBeagleDataLikelihoodDelegate(benchmark.likelihood);
        }
    }
}
