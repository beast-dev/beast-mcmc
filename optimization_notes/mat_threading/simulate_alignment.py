#!/usr/bin/env python3
"""
Makes a MAT benchmark XML with an alignment simulated on a given tree, for the profiling study of the
MAT's cost on large alignments (many variable sites, many events per branch).

The template is a benchmark MAT XML (e.g. benchmark-big-3-mat.xml): its alignment is replaced by sequences
simulated under HKY (kappa 2, equal frequencies) with a strict clock at RATE substitutions per site per
year along TREE (newick, branch lengths in years; the tree the XML's coalescent simulator produced), and
the coalescent starting tree by that tree, so that the chain starts at the truth. The variable fraction of
sites is printed: it is what the MAT's cost depends on, not the number of sites.

usage: simulate_alignment.py TEMPLATE_XML TREE_NEWICK OUT_XML --sites N --rate R [--kappa 2] [--seed 1]
       [--chain-length N] [--log-every N] [--name NAME]
"""
import argparse
import re
import sys

import numpy as np


def parse_newick(text):
    """@return (children[node], length[node], label[node], root): a tree as parallel lists."""
    text = re.sub(r'\[&[^\]]*\]', '', text.strip())
    children, lengths, labels = [], [], []
    stack = []
    i = 0
    n = len(text)

    def new_node():
        children.append([])
        lengths.append(0.0)
        labels.append(None)
        return len(children) - 1

    root = None
    current = None
    while i < n:
        c = text[i]
        if c == '(':
            node = new_node()
            if stack:
                children[stack[-1]].append(node)
            stack.append(node)
            current = None
            i += 1
        elif c == ',':
            current = None
            i += 1
        elif c == ')':
            current = stack.pop()
            i += 1
        elif c == ':':
            m = re.match(r':([0-9.eE+-]+)', text[i:])
            lengths[current] = float(m.group(1))
            i += len(m.group(0))
        elif c == ';':
            root = current
            break
        elif c in ' \t\n\r':
            i += 1
        else:
            if c in "'\"":
                # a quoted label (BEAST double-quotes labels with spaces, without escaping inner quotes:
                # the label ends at the quote before the branch length)
                j = text.index(c + ':', i + 1)
                label = text[i + 1:j]
                i = j + 1
            else:
                m = re.match(r"[^,():;]+", text[i:])
                label = m.group(0)
                i += len(label)
            if current is None:
                current = new_node()
                children[stack[-1]].append(current)
            labels[current] = label
    return children, lengths, labels, root


def to_newick(children, lengths, labels, root):
    """@return the tree as newick with every label single-quoted (inner quotes doubled)"""
    parts = []

    def write(node):
        if children[node]:
            parts.append('(')
            for i, child in enumerate(children[node]):
                if i > 0:
                    parts.append(',')
                write(child)
            parts.append(')')
        else:
            parts.append("'" + labels[node].replace("'", "''") + "'")
        if node != root:
            parts.append(':%.17g' % lengths[node])

    import sys
    sys.setrecursionlimit(max(10000, sys.getrecursionlimit()))
    write(root)
    return ''.join(parts) + ';'


def hky_matrix(kappa, freqs):
    """@return Q for HKY, scaled to one expected substitution per site per unit distance"""
    s = len(freqs)
    q = np.zeros((s, s))
    for i in range(s):
        for j in range(s):
            if i != j:
                transition = (i, j) in ((0, 2), (2, 0), (1, 3), (3, 1))  # A<->G, C<->T
                q[i, j] = freqs[j] * (kappa if transition else 1.0)
        q[i, i] = -q[i].sum()
    scale = -sum(freqs[i] * q[i, i] for i in range(s))
    return q / scale


def transition_probabilities(q, distance):
    w, v = np.linalg.eig(q)
    p = (v @ np.diag(np.exp(w * distance)) @ np.linalg.inv(v)).real
    p = np.clip(p, 0.0, None)
    return p / p.sum(axis=1, keepdims=True)


def simulate(children, lengths, labels, root, sites, rate, kappa, rng):
    freqs = np.array([0.25, 0.25, 0.25, 0.25])
    q = hky_matrix(kappa, freqs)
    sequences = {}
    root_states = rng.choice(4, size=sites, p=freqs).astype(np.int8)
    stack = [(root, root_states)]
    while stack:
        node, states = stack.pop()
        for child in children[node]:
            p = transition_probabilities(q, lengths[child] * rate)
            child_states = np.empty(sites, dtype=np.int8)
            for a in range(4):
                idx = np.nonzero(states == a)[0]
                if len(idx) > 0:
                    child_states[idx] = rng.choice(4, size=len(idx), p=p[a])
            if children[child]:
                stack.append((child, child_states))
            else:
                sequences[labels[child]] = child_states
    return sequences


def xml_escape(label):
    return label.replace('&', '&amp;').replace('<', '&lt;').replace('"', '&quot;')


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('expected exactly one occurrence of: ' + old[:80])
    return text.replace(old, new)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('template')
    parser.add_argument('tree')
    parser.add_argument('out')
    parser.add_argument('--sites', type=int, required=True)
    parser.add_argument('--rate', type=float, required=True)
    parser.add_argument('--kappa', type=float, default=2.0)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--chain-length', type=int, default=1000000)
    parser.add_argument('--log-every', type=int, default=10000)
    parser.add_argument('--name', default=None)
    args = parser.parse_args()
    name = args.name or re.sub(r'\.xml$', '', args.out.split('/')[-1])

    with open(args.tree) as f:
        newick = f.read()
    children, lengths, labels, root = parse_newick(newick)
    tips = [l for l in labels if l is not None]
    with open(args.template) as f:
        template_ids = set(re.findall(r'<taxon id="([^"]*)"', f.read()))
    missing = [t for t in tips if xml_escape(t) not in template_ids]
    if missing:
        raise ValueError('%d tree labels are not taxon ids of the template, e.g. %s' % (len(missing), missing[:3]))
    rng = np.random.default_rng(args.seed)
    sequences = simulate(children, lengths, labels, root, args.sites, args.rate, args.kappa, rng)
    matrix = np.stack([sequences[t] for t in tips])
    variable = int(np.count_nonzero((matrix != matrix[0]).any(axis=0)))
    print('%s: %d taxa, %d sites, %d variable (%.1f%%)' % (name, len(tips), args.sites, variable,
                                                            100.0 * variable / args.sites), file=sys.stderr)

    with open(args.template) as f:
        xml = f.read()
    # the alignment
    start = xml.index('\t<alignment id="alignment" dataType="nucleotide">')
    end = xml.index('\t</alignment>', start) + len('\t</alignment>')
    codes = np.array(list('ACGT'))
    rows = ['\t<alignment id="alignment" dataType="nucleotide">']
    for t in tips:
        rows.append('\t\t<sequence>\n\t\t\t<taxon idref="%s"/>\n\t\t\t%s\n\t\t</sequence>' %
                    (xml_escape(t), ''.join(codes[sequences[t]])))
    rows.append('\t</alignment>')
    xml = xml[:start] + '\n'.join(rows) + xml[end:]
    # the starting tree: the simulation's tree
    start = xml.index('\t<coalescentSimulator id="startingTree">')
    end = xml.index('\t</coalescentSimulator>', start) + len('\t</coalescentSimulator>')
    newick_text = to_newick(children, lengths, labels, root).replace('&', '&amp;').replace('<', '&lt;')
    xml = xml[:start] + '\t<newick id="startingTree" units="years">\n\t\t%s\n\t</newick>' % newick_text + xml[end:]
    xml = replace_once(xml, '<coalescentTree idref="startingTree"/>', '<newick idref="startingTree"/>')
    # no ambiguities: plain constant-site pooling
    xml = re.sub(r'(<patterns id="patterns"[^>]*?) ambiguousConstant="true" ambiguityThreshold="[0-9.]+"', r'\1', xml)
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % args.chain_length, xml)
    xml = re.sub(r'(<log id="(?:screenLog|fileLog)" logEvery=)"\d+"', r'\1"%d"' % args.log_every, xml)
    xml = re.sub(r'(<logTree id="treeFileLog" logEvery=)"\d+"', r'\1"%d"' % max(args.log_every, args.chain_length // 10), xml)
    xml = re.sub(r'fileName="[\w-]+\.log"', 'fileName="%s.log"' % name, xml)
    xml = re.sub(r'fileName="[\w-]+\.trees"', 'fileName="%s.trees"' % name, xml)
    xml = re.sub(r'<!-->.*?<-->', '<!-->\nSimulated: %d sites under HKY (kappa %g) at rate %g on the big-3 coalescent tree; '
                 '%d variable sites (%.1f%%); seed %d\n<-->' % (args.sites, args.kappa, args.rate, variable,
                                                              100.0 * variable / args.sites, args.seed), xml,
                 count=1, flags=re.S)
    with open(args.out, 'w') as f:
        f.write(xml)
    print('wrote', args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
