#!/usr/bin/env python3
"""
Fitch parsimony score of trees from a BEAST .trees file, for the alignment in a BEAST XML.
Ambiguous states (anything but A, C, G, T) are treated as missing, as the MAT does.

usage: tree_parsimony.py XML TREES [every]
"""
import html
import re
import sys

import numpy as np

CODE = {'A': 1, 'C': 2, 'G': 4, 'T': 8}


def read_alignment(xml):
    text = open(xml).read()
    seqs = {}
    for m in re.finditer(r'<sequence>\s*<taxon idref="([^"]+)"/>\s*([^<]+)</sequence>', text):
        seqs[m.group(1)] = re.sub(r'\s+', '', m.group(2)).upper()
    return seqs


def read_trees(path):
    translate = {}
    trees = []
    in_translate = False
    for line in open(path):
        s = line.strip()
        if s.lower() == 'translate':
            in_translate = True
            continue
        if in_translate:
            if s == ';':
                in_translate = False
                continue
            number, name = s.rstrip(',').split(None, 1)
            translate[number] = name.strip("'\"")
        elif s.lower().startswith('tree '):
            state = int(re.search(r'STATE_(\d+)', s).group(1))
            trees.append((state, s[s.index('('):].rstrip(';')))
    return translate, trees


def parse_newick(newick):
    """@return children lists and tip names, as nested lists (iterative)."""
    newick = re.sub(r'\[[^\]]*\]', '', newick)
    newick = re.sub(r':[-\d.eE+]+', '', newick)
    stack = [[]]
    token = ''
    for ch in newick:
        if ch == '(':
            stack.append([])
        elif ch in ',)':
            if token:
                stack[-1].append(token)
                token = ''
            if ch == ')':
                node = stack.pop()
                stack[-1].append(node)
        else:
            token += ch
    if token:
        stack[-1].append(token)
    return stack[0][0]


def fitch(tree, sets_by_tip):
    """@return the Fitch score of the (nested-list) tree, summed over sites (vectorised)."""
    score = 0

    def visit(node):
        nonlocal score
        if isinstance(node, str):
            return sets_by_tip[node]
        child_sets = [visit(c) for c in node]
        current = child_sets[0]
        for other in child_sets[1:]:
            # empty sets (missing data) are ignored
            both = (current != 0) & (other != 0)
            inter = current & other
            union = current | other
            score += int(np.count_nonzero(both & (inter == 0)))
            current = np.where(both, np.where(inter != 0, inter, union), current | other)
        return current

    sys.setrecursionlimit(100000)
    visit(tree)
    return score


def read_matrix(xml):
    """@return the taxon names and their state sets (bit codes; 0 for missing) at the variable sites."""
    seqs = read_alignment(xml)
    names = list(seqs)
    length = len(seqs[names[0]])
    matrix = np.zeros((len(names), length), dtype=np.uint8)
    for i, name in enumerate(names):
        matrix[i] = [CODE.get(c, 0) for c in seqs[name]]
    # only sites with at least two different observed states matter
    variable = np.array([len(set(matrix[:, j][matrix[:, j] != 0])) > 1 for j in range(length)])
    return names, matrix[:, variable]


def main():
    xml, trees_path = sys.argv[1], sys.argv[2]
    every = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    names, matrix = read_matrix(xml)
    sets_by_name = {key(name): matrix[i] for i, name in enumerate(names)}
    translate, trees = read_trees(trees_path)
    for state, newick in trees[::every]:
        tree = parse_newick(newick)
        tree = relabel(tree, translate)
        print(state, fitch(tree, sets_by_name))


def key(name):
    return re.sub(r'[\'"\s]', '', html.unescape(name))


def relabel(node, translate):
    if isinstance(node, str):
        return key(translate.get(node, node))
    return [relabel(c, translate) for c in node]


if __name__ == '__main__':
    main()
