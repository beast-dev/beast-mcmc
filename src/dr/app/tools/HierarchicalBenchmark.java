/*
 * Copyright (c) 2026 the BEAST Development Team.
 * Licensed under the GNU Lesser General Public License, version 2 or later.
 */
package dr.app.tools;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.InstanceDetails;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.MessageDigest;
import java.util.*;

/** Direct native BEAGLE driver for Jeff's four paired likelihood calculations. */
public final class HierarchicalBenchmark {
    private enum Method { H, HA, NE, AE }

    private static final class Options {
        Path input, output;
        Method method;
        String endpoint = "pruning";
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
            if (!o.endpoint.equals("pruning") && !o.endpoint.equals("model-update"))
                throw new IllegalArgumentException("--endpoint must be pruning or model-update");
            if (o.warmupEvals < 0 || o.warmupSeconds < 0 || o.batchSeconds <= 0 || o.batches < 1
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

    private static final class Phases {
        long model, matrices, upload, pruning, root, evaluations, eigensystems, modelUploads, matrixUpdates;
        void clear() { model = matrices = upload = pruning = root = evaluations = eigensystems = modelUploads = matrixUpdates = 0; }
        Map<String, Object> json() {
            Map<String, Object> p = map();
            p.put("evaluations", evaluations); p.put("eigensystem_constructions", eigensystems);
            p.put("model_uploads", modelUploads);
            p.put("updated_transition_matrices", matrixUpdates);
            p.put("model_seconds", model / 1e9); p.put("matrix_seconds", matrices / 1e9);
            p.put("data_upload_seconds", upload / 1e9); p.put("pruning_seconds", pruning / 1e9);
            p.put("root_seconds", root / 1e9);
            return p;
        }
    }

    private final Options options;
    private final Input input;
    private final HierarchicalSubstitutionModel model;
    private final Parameter rates;
    private final int rootRateIndex, blockSize, blockCount, padding;
    private final boolean streaming, hierarchical;
    private final Beagle beagle;
    private final int[][] blockTips;
    private final double[] blockWeights, rootResult = new double[1], siteBuffer;
    private final int[] rootIndex, zero = {0}, none = {Beagle.NONE};
    private final Phases phases = new Phases();
    private int point;
    private double[] collected;

    private HierarchicalBenchmark(Options options, Input input) throws Exception {
        this.options = options; this.input = input;
        blockSize = options.blockSize == 0 ? input.k : Math.min(options.blockSize, input.k);
        blockCount = (input.k + blockSize - 1) / blockSize;
        padding = blockCount * blockSize - input.k;
        streaming = blockCount > 1;
        hierarchical = options.method == Method.H || options.method == Method.HA;
        String[] labels = new String[input.s];
        for (int i = 0; i < input.s; i++) labels[i] = String.format(Locale.ROOT, "state%03d", i);
        TreeModel tree = new DefaultTreeModel("characterTree", new NewickImporter(input.newick(input.root) + ";").importTree(null));
        rates = new Parameter.Default(input.s - 1, 1.0);
        fillRates(tree, tree.getRoot(), input.ratesByClade());
        rootRateIndex = tree.getRoot().getNumber() - input.s;
        FrequencyModel frequency = new FrequencyModel(new GeneralDataType(labels), new Parameter.Default(input.pi.clone()));
        model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, tree, frequency, true));
        model.setEigenSystemType(options.method == Method.NE
                ? HierarchicalSubstitutionModel.EigenSystemType.NUMERICAL : HierarchicalSubstitutionModel.EigenSystemType.ANALYTICAL);
        long flags = BeagleFlag.PROCESSOR_CPU.getMask() | BeagleFlag.PRECISION_DOUBLE.getMask()
                | BeagleFlag.THREADING_NONE.getMask() | BeagleFlag.SCALING_MANUAL.getMask();
        if (hierarchical) flags |= BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()
                | (options.method == Method.HA ? BeagleFlag.HIERARCHICAL_AGGREGATE_ON.getMask()
                : BeagleFlag.HIERARCHICAL_AGGREGATE_OFF.getMask());
        else flags |= BeagleFlag.EIGEN_REAL.getMask() | BeagleFlag.VECTOR_SSE.getMask();
        beagle = BeagleFactory.loadBeagleInstance(input.n, 2 * input.n - 1, input.n, input.s,
                blockSize, 1, input.matrices.length, 1, 0, null, flags, flags);
        InstanceDetails details = beagle.getDetails();
        if (details == null || (details.getFlags() & flags) != flags)
            throw new IllegalStateException("Required native implementation flags unavailable: " + BeagleFlag.toString(flags));
        if (!hierarchical && (details.getFlags() & (BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()
                | BeagleFlag.SPECTRAL_REPRESENTATION.getMask())) != 0)
            throw new IllegalStateException("Dense arm selected a nonstandard representation");
        beagle.setStateFrequencies(0, input.pi); beagle.setCategoryRates(new double[]{1.0});
        beagle.setCategoryWeights(0, new double[]{1.0});
        rootIndex = new int[]{2 * input.n - 2};
        blockTips = new int[input.n][blockSize]; blockWeights = new double[blockSize]; siteBuffer = new double[blockSize];
        uploadBlock(0);
        uploadModel();
    }

    private BitSet fillRates(TreeModel tree, NodeRef node, Map<BitSet, Double> byClade) {
        BitSet clade = new BitSet(input.s);
        if (tree.isExternal(node)) clade.set(Integer.parseInt(tree.getTaxonId(node.getNumber()).substring(5)));
        else {
            for (int c = 0; c < tree.getChildCount(node); c++) clade.or(fillRates(tree, tree.getChild(node, c), byClade));
            Double rate = byClade.get(clade);
            if (rate == null) throw new IllegalArgumentException("Character tree/rate mapping mismatch: " + clade);
            rates.setParameterValue(node.getNumber() - input.s, rate);
        }
        return clade;
    }

    private void uploadBlock(int offset) {
        long start = System.nanoTime();
        int actual = Math.min(blockSize, input.k - offset);
        for (int tip = 0; tip < input.n; tip++) {
            System.arraycopy(input.tips[tip], offset, blockTips[tip], 0, actual);
            // Duplicate this block's first real pattern; zero multiplicity makes padding likelihood-neutral.
            Arrays.fill(blockTips[tip], actual, blockSize, blockTips[tip][0]);
            beagle.setTipStates(tip, blockTips[tip]);
        }
        System.arraycopy(input.weights, offset, blockWeights, 0, actual);
        Arrays.fill(blockWeights, actual, blockSize, 0.0);
        beagle.setPatternWeights(blockWeights);
        phases.upload += System.nanoTime() - start;
    }

    private void uploadModel() {
        long start = System.nanoTime();
        if (hierarchical) beagle.setEigenDecomposition(0, model.getBeagleModelTriples(),
                model.getBeagleModelUnused(), model.getBeagleModelLambdas());
        else {
            EigenDecomposition eigen = model.getEigenDecomposition();
            beagle.setEigenDecomposition(0, eigen.getEigenVectors(), eigen.getInverseEigenVectors(), eigen.getEigenValues());
            phases.eigensystems++;
        }
        phases.modelUploads++;
        phases.model += System.nanoTime() - start;
        start = System.nanoTime();
        beagle.updateTransitionMatrices(0, input.matrices, null, null, input.matrixLengths, input.matrices.length);
        phases.matrixUpdates += input.matrices.length;
        phases.matrices += System.nanoTime() - start;
    }

    private void setPoint(int next) {
        if (next != point) {
            long start = System.nanoTime();
            rates.setParameterValue(rootRateIndex, input.rawRates[input.root] * (next == 0 ? 1.0 : 1.01));
            phases.model += System.nanoTime() - start;
            point = next;
            uploadModel();
        }
    }

    private double evaluate(boolean update, boolean collect) {
        if (update) setPoint(1 - point);
        if (collect) collected = new double[input.k];
        double total = 0.0;
        for (int block = 0; block < blockCount; block++) {
            int offset = block * blockSize;
            if (streaming) uploadBlock(offset);
            long start = System.nanoTime();
            beagle.updatePartials(input.operations, input.n - 1, Beagle.NONE);
            phases.pruning += System.nanoTime() - start;
            start = System.nanoTime();
            beagle.calculateRootLogLikelihoods(rootIndex, zero, zero, none, 1, rootResult);
            phases.root += System.nanoTime() - start;
            if (!Double.isFinite(rootResult[0])) throw new ArithmeticException("Nonfinite native likelihood");
            total += rootResult[0];
            if (collect) {
                beagle.getSiteLogLikelihoods(siteBuffer);
                System.arraycopy(siteBuffer, 0, collected, offset, Math.min(blockSize, input.k - offset));
            }
        }
        phases.evaluations++;
        return total;
    }

    private static void saveDoubles(Path path, double[] values) throws IOException {
        try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(path)))) {
            for (double value : values) {
                if (!Double.isFinite(value)) throw new ArithmeticException("Nonfinite per-pattern log likelihood");
                out.writeDouble(value);
            }
        }
    }

    private Map<String, Object> run(long launchStart) throws IOException {
        Map<String, Object> result = map();
        double base = evaluate(false, true);
        double preparationSeconds = (System.nanoTime() - launchStart) / 1e9;
        Map<String, Object> first = map(); first.put("log_likelihood", base);
        first.put("java_main_to_first_seconds", preparationSeconds);
        System.out.println("FIRST_LIKELIHOOD " + json(first)); System.out.flush();
        Path baseFile = options.output.resolveSibling(options.output.getFileName() + ".base.f64be");
        Path changedFile = options.output.resolveSibling(options.output.getFileName() + ".changed.f64be");
        saveDoubles(baseFile, collected);
        setPoint(1);
        double changed = evaluate(false, true); saveDoubles(changedFile, collected);
        setPoint(0); evaluate(false, false);
        Map<String, Object> validation = map();
        validation.put("base_log_likelihood", base); validation.put("changed_log_likelihood", changed);
        validation.put("base_site_logs", baseFile.toString()); validation.put("changed_site_logs", changedFile.toString());
        validation.put("site_log_dtype", ">f8");
        result.put("validation", validation);
        result.put("preparation_seconds", preparationSeconds);
        result.put("preparation_clock_origin", "Java main entry; timestamp FIRST_LIKELIHOOD stdout marker for process-launch timing");
        result.put("first_log_likelihood", base);
        result.put("method", options.method.name()); result.put("endpoint", options.endpoint);
        result.put("states", input.s); result.put("tips", input.n); result.put("sites", input.l); result.put("patterns", input.k);
        result.put("execution_mode", streaming ? "streamed" : "resident");
        result.put("native_patterns", blockSize); result.put("pattern_blocks", blockCount); result.put("padding_patterns", padding);
        result.put("shared_transition_matrices", input.sharedTransitions);
        result.put("matrix_buffer_count", input.matrices.length);
        result.put("distinct_branch_lengths", input.distinctBranchLengths);
        result.put("transition_reuse_note", "Each phylogeny branch has its own transition matrix, including when branch lengths are equal.");
        result.put("scaling", "manual, no rescaling operations; all finite likelihoods required");
        result.put("category_count", 1); result.put("threads", 1);
        result.put("root_rate_multiplier_base", 1.0); result.put("root_rate_multiplier_changed", 1.01);
        result.put("input", options.input.toString()); result.put("input_sha256", sha256(options.input));
        InstanceDetails details = beagle.getDetails();
        String implementation = details.getImplementationName();
        result.put("implementation", implementation == null ? (hierarchical ? "CPU-Hierarchical-Double"
                : input.s == 4 ? "CPU-4State-SSE-Double" : "CPU-SSE-Double") : implementation);
        result.put("implementation_name_source", implementation == null ? "inferred from strictly required native flags" : "native JNI");
        result.put("flags", details.getFlags()); result.put("flag_names", BeagleFlag.toString(details.getFlags()));
        result.put("resource", details.getResourceNumber());
        result.put("aggregation_requested", options.method == Method.HA ? "on" : options.method == Method.H ? "off" : "not applicable");
        result.put("aggregation_effective_note", options.method == Method.HA && input.s == 8
                ? "current S8 shape dispatch retains full state kernel" : "see selected native flags and verbose log");
        result.put("streaming_note", streaming
                ? "Every complete evaluation includes tip/weight uploads and any aggregation class rebuilds for all blocks. Model/eigensystem and transition updates occur once per Q, not per block."
                : "Patterns and transition data remain resident for cached-pruning evaluations.");
        List<Object> batches = new ArrayList<>();
        result.put("batches", batches);
        if (!options.validateOnly) {
            boolean updates = options.endpoint.equals("model-update");
            long start = System.nanoTime(), count = 0;
            do {
                evaluate(updates, false); count++;
                if (count >= options.maxEvals && (count < options.warmupEvals
                        || (System.nanoTime() - start) / 1e9 < options.warmupSeconds))
                    throw new IllegalStateException("Warmup max-evals reached before requested duration");
            } while (count < options.warmupEvals || (System.nanoTime() - start) / 1e9 < options.warmupSeconds);
            Map<String, Object> warmup = map(); warmup.put("evaluations", count);
            warmup.put("seconds", (System.nanoTime() - start) / 1e9); result.put("warmup", warmup);
            double[] times = new double[options.batches];
            for (int b = 0; b < options.batches; b++) {
                phases.clear(); start = System.nanoTime(); count = 0;
                double last;
                do {
                    last = evaluate(updates, false); count++;
                    if (count >= options.maxEvals && (System.nanoTime() - start) / 1e9 < options.batchSeconds)
                        throw new IllegalStateException("Batch max-evals reached before requested duration");
                } while ((System.nanoTime() - start) / 1e9 < options.batchSeconds || (updates && (count & 1) != 0));
                double seconds = (System.nanoTime() - start) / 1e9;
                Map<String, Object> batch = map();
                batch.put("batch", b + 1); batch.put("iterations", count); batch.put("seconds", seconds);
                batch.put("ms_per_eval", seconds * 1000.0 / count); batch.put("last_log_likelihood", last);
                batch.put("phases", phases.json()); batch.put("phylogeny_operations", count * (input.n - 1) * blockCount);
                batches.add(batch); times[b] = seconds * 1000.0 / count;
            }
            Arrays.sort(times);
            result.put("median_ms_per_eval", times.length % 2 == 1 ? times[times.length / 2]
                    : 0.5 * (times[times.length / 2 - 1] + times[times.length / 2]));
        }
        result.put("status", "ok"); result.put("format_version", 1);
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
        long launchStart = System.nanoTime(); Options options = null; HierarchicalBenchmark benchmark = null;
        try {
            options = Options.parse(args); Files.createDirectories(options.output.getParent());
            if (Boolean.getBoolean("java.only")) throw new IllegalArgumentException("Native BEAGLE is required; java.only is not supported");
            if (System.getenv("BEAGLE_HRATE_AGGREGATE") != null)
                throw new IllegalArgumentException("Unset BEAGLE_HRATE_AGGREGATE: method selection uses native flags");
            Input input = Input.read(options.input);
            benchmark = new HierarchicalBenchmark(options, input);
            Map<String, Object> result = benchmark.run(launchStart); saveJson(options.output, result);
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
            if (benchmark != null) benchmark.beagle.finalize();
        }
    }
}
