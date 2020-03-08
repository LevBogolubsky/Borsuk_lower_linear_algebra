#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import sys

DEFAULT_PRIME_POWERS_FILENAME = 'prime_powers.tsv'

TEX_BEGIN = """
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}

\usepackage{geometry}
 \geometry{
 a4paper,
 total={170mm,257mm},
 left=5mm,
 top=20mm,
 }

\usepackage{longtable}
\\begin{document}
\\begin{center}
\\scriptsize
"""

TEX_END = """
\\normalsize
\hspace*{-1cm}
\end{center}
\end{document}
"""


def load_prime_powers(filename):
    return [int(line.strip('\n')) for line in open(filename)]


def print_table(dim_from, dim_to, systems, systems_to_params, result, good_dims, fout):
    def format_system_header_part(s):
        return ''.join(['\\multicolumn{2} {|c|} {', s, '}'])

    print('\\begin{longtable}{|c|' + 'c|c|' * len(systems) * 2 + '}', file=fout)
    print('\\caption{Таблица \\label{Table}} \\\\', file=fout)

    print('\\hline', file=fout)
    print('серия графов $\\rightarrow$ & ' + ' & '.join(map(format_system_header_part, systems)) + ' \\\\ \hline',  file=fout)

    parameter_header_parts = []
    for system in systems:
        parameter_header_parts.extend(['$f$', systems_to_params[system]])

    print('размерность $\\downarrow$ & ' + ' & '.join(parameter_header_parts) + ' \\\\ \hline',  file=fout)
    print('\\endhead', file=fout)

    for dim in range(dim_from, dim_to + 1):
        results = map(lambda x : str(result[dim][x]), systems)
        if dim in good_dims:
            print(dim, ' & ' + ' & '.join(results) + ' \\\\ \hline',  file=fout)

    print('\\end{longtable}', file=fout)


def make_bold(dim_from, dim_to, systems, result):
    for dim in range(dim_from, dim_to + 1):
        max_value = 1

        for system in systems:
            max_value = max(max_value, result[dim][system].estimate)

        for system in systems:
            if result[dim][system].estimate == max_value:
                result[dim][system].is_max = True


def update_progress_bar(count, total):
    PROGRESS_BAR_SIZE = 100 
    filled_part_size = int(round(PROGRESS_BAR_SIZE * count / float(total)))
    
    bar = '|' * filled_part_size + '-' * (PROGRESS_BAR_SIZE - filled_part_size)
    print(
        '[{0}] {1}{2}\r'.format(
            bar, 
            int(round(100.0 * count / float(total), 1)),
            '%'
        ),
        file=sys.stderr,
        end=''
    )
    sys.stderr.flush()

    if count == total:
        print('\n', file=sys.stderr, end='')
        sys.stderr.flush()
