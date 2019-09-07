#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from utils import (
    DEFAULT_PRIME_POWERS_FILENAME,
    TEX_BEGIN,
    TEX_END,
    load_prime_powers,
    print_table,
    make_bold,
    update_progress_bar,
)

from argparse import ArgumentParser
from math import ceil
from collections import defaultdict
from scipy.special import binom
from primefac import primefac

import sys 

DEFAULT_LOWEST_COORDINATE_COUNT = 5
DEFAULT_HIGHEST_COORDINATE_COUNT = 45

DEFAULT_LOWEST_DIMENSION = 1
DEFAULT_HIGHEST_DIMENSION = 1000

SYSTEM_G1_THROW_OUT ='$G_1$'
SYSTEM_G2_THROW_OUT = '$G_2$'

SYSTEM_MOD_0_FREE = '$G_3$'
SYSTEM_MOD_MINUS_1_FREE = '$G_4$'
SYSTEM_MOD_PLUS_1_FREE = '$G_5$'

SYSTEM_TWO_RESTRICTIONS = '$G_6$'
SYSTEM_TWO_RESTRICTIONS_MINUS_1 = '$G_7$'
SYSTEM_TWO_RESTRICTIONS_MINUS_2 = '$G_8$'
SYSTEM_TWO_RESTRICTIONS_MINUS_3 = '$G_9$'


class BorsukLowerEstimate(object):
    def __init__(self, outer_estimate=None, outer_parameters=None):
        if outer_estimate is None:
            self.estimate = 1
        else:
            self.estimate = outer_estimate
        if outer_parameters is None:
            self.parameters = None
        else:
            self.parameters = outer_parameters
        self.is_max = False

    def __str__(self):
        result = None
        if self.parameters is not None:
            if len(self.parameters) == 1: 
                parameters_string = '({0})'.format(self.parameters[0])
            else:
                parameters_string = str(self.parameters)
            result = '\t'.join([str(self.estimate), parameters_string])
        else:
            result = str(self.estimate)
        if self.is_max:
            result = '\\textbf{' + result + '}'
        return result

    def is_from_spherical_enrichment(self):
        return self.parameters is None

    def update(self, outer_estimate, outer_parameters):
        self.estimate = outer_estimate
        self.parameters = outer_parameters

    def update_if_better(self, outer_estimate):
        if outer_estimate.estimate > self.estimate:
            self.estimate = outer_estimate.estimate
            self.parameters = outer_estimate.parameters


def try_update_affected_dimensions(result, dimension, graph_type, affected_dimensions):
    if result[dimension][graph_type].estimate > dimension + 1:
        affected_dimensions.add(dimension)
    return affected_dimensions


# Lemma 16
def calculate_dimension_simple_dual_transition(n): 
    return n * (n + 1) / 2 - 1


# Lemma 15, claim 10
def calculate_dimension_two_digits_single_ban_custom_remainder(n): 
    return int(binom(n, 2))


# Lemma 15, claims 10 and 13
def calculate_dimension_two_restrictions(n, I, J, F):
    assert F in [J, J - 1]
    return int(binom(n - J + 1, 2) + (n - J + 1) * (J - F))


# Theorem 1
def calculate_alpha_upper_estimate(n, q, digits):
    result = 0

    # Theorem 4
    if digits == 3:
        for m_1 in range(n + 1):
            for m_2 in range(n + 1):
                if m_1 + 2 * m_2 in [q - 1, q - 2]: 
                    result += binom(n, m_1) * binom(n - m_1, m_2)
    elif digits == 2:
        result = binom(n, q - 1)
    
    return int(result)


# Subsection 6.1
def find_f_best_lower_estimate_fixed(n, prime_powers):
    best_answer = BorsukLowerEstimate()
    
    for q in prime_powers:
        if q - 1 > n:
            break
        if q % 2 == 0:
            continue
        alpha_upper = calculate_alpha_upper_estimate(n, q, 3)
        d = q 
        if d < n:
            k_0 = n - d
            for k_1 in range(1, d):
                k_m1 = d - k_1
                size = int(binom(n, k_0) * binom(n - k_0 - 1, k_m1))
                f_estimate = int(ceil(float(size) / alpha_upper))
                if best_answer.estimate <= f_estimate:
                    best_answer.update(f_estimate, (k_m1, k_0, k_1))
    
    return best_answer


# Subsection 6.2
def find_f_best_lower_estimate_free(n, prime_powers): 
    best_answer = BorsukLowerEstimate()
    
    for q in prime_powers:
        if q - 1 > n:
            break
        if q % 2 == 0:
            continue
        alpha_upper = calculate_alpha_upper_estimate(n, q, 3)
        d = q 
        if d < n:
            size = int(binom(n, d) * (2 ** (d - 1))) 
            f_estimate = int(ceil(float(size) / alpha_upper))
            if best_answer.estimate <= f_estimate:
                best_answer.update(f_estimate, (n, d))

    return best_answer


def get_remainder_modulo_four_group(n):
    """
    group 0 -> epsilon = 0, -4, ...
    group 1 -> epsilon = 1, -3, ...
    group 2 -> epsilon = 2, -2, ...
    group 3 -> epsilon = -1, -5, ...
    """
    group = 4 - (n % 4) if n % 4 != 0 else 0 
    assert group in [0, 1, 2, 3]
    return group


# Subsection 6.3
def find_f_two_digits_single_ban_custom_remainder(n, prime_powers):
    best_answer = BorsukLowerEstimate()

    group = get_remainder_modulo_four_group(n)
    
    epsilon = None
    if group == 0:
        epsilon = 0
    elif group == 1:
        epsilon = 1
    elif group == 3:
        epsilon = -1
    else:
        raise Exception

    q = (n + epsilon) / 4

    if q not in prime_powers:
        return best_answer
    if q % 2 == 0:
        return best_answer
    
    if epsilon in [0, 1]:
        size = 2 ** (n - 2)
    else:
        assert epsilon == -1
        size = 2 ** (n - 3)

    if epsilon in [0, 1]:
        alpha_upper = binom(n - 1, q - 1)
    else:
        alpha_upper = binom(n - 2, q - 1)

    f_estimate = int(ceil(float(size) / alpha_upper))

    if best_answer.estimate <= f_estimate:
        best_answer.update(f_estimate, (n,))

    return best_answer    


# Subsection 6.4
def find_f_two_restrictions_common_remainder(n, prime_powers):
    RESTRICTION_MAPPING = {
        0: (0, 4, 3),
        -4: (4, 8, 7),
        1: (-1, 3, 2),
        -3: (3, 7, 6),
        2: (-2, 2, 1),
        -2: (2, 6, 5),
        -1: (1, 5, 4),
        -5: (5, 9, 8),
    }

    SYSTEM_MAPPING = {  
        0: SYSTEM_TWO_RESTRICTIONS,
        -4: SYSTEM_TWO_RESTRICTIONS,
        1: SYSTEM_TWO_RESTRICTIONS_MINUS_1,
        -3: SYSTEM_TWO_RESTRICTIONS_MINUS_1,
        2: SYSTEM_TWO_RESTRICTIONS_MINUS_2,
        -2: SYSTEM_TWO_RESTRICTIONS_MINUS_2,
        -1: SYSTEM_TWO_RESTRICTIONS_MINUS_3,
        -5: SYSTEM_TWO_RESTRICTIONS_MINUS_3,
    }

    K = 1

    group = get_remainder_modulo_four_group(n)

    q = None
    epsilon = None

    if group == 0 and n % 4 == 0:
        if n / 4 in prime_powers:
            q = n / 4
            epsilon = 0
        else:
            q = (n - 4)  / 4 
            epsilon = -4 
    elif group == 1 and (n + 1) % 4 == 0: 
        if (n + 1) / 4 in prime_powers:
            q = (n + 1) / 4
            epsilon = 1
        else:
            q = (n - 3) / 4
            epsilon = -3
    elif group == 2 and (n + 2) % 4 == 0:
        if (n + 2) / 4 in prime_powers:
            q = (n + 2) / 4
            epsilon = 2
        else:
            q = (n - 2) / 3
            epsilon = -2
    elif group == 3 and (n - 1) % 4 == 0:
        if (n - 1) / 4 in prime_powers:
            q = (n - 1) / 4
            epsilon = -1
        else:
            q = (n - 5) / 4
            epsilon = -5
    else:
        raise Exception

    if q not in prime_powers:
        return
    if q % 2 == 0:
        return
    if K >= q - 1:
        return

    r_1, r_2, min_F = RESTRICTION_MAPPING[epsilon]  
    system = SYSTEM_MAPPING[epsilon]

    if (r_1 - r_2) % primefac(q).next() == 0:
        return

    I = r_1 + 1 
    J = r_2 + 1

    while I <= J:
        if min(I, J) >= 1:
            for F in range(I, J + 1):
                if (F >= I) and (F >= J - 1) and (F <= n - 1) and (F >= min_F) and ((F == J) or (F == J - 1)):
                    dimension = calculate_dimension_two_restrictions(n, I, J, F)
            
                    size = 2 ** (n - F - 1)
                    alpha_upper = 0

                    for k in range(q - K):
                        alpha_upper += binom(n - F, k)

                    f_min = int(ceil(float(size) / alpha_upper))

                    if f_min > dimension + 1:
                        yield system, dimension, BorsukLowerEstimate(f_min, (n, I, J, F))
        I += 1
        J -= 1


def perform_spherical_enrichment(
        lowest_dimension, 
        highest_dimension, 
        systems, 
        result, 
        affected_dimensions
):
    for dimension in range(lowest_dimension, highest_dimension + 1):
        max_value = 1

        for system in systems:
            if result[dimension][system].estimate < result[dimension - 1][system].estimate + 1:
                result[dimension][system] = BorsukLowerEstimate(result[dimension - 1][system].estimate + 1, None)
            max_value = max(max_value, result[dimension][system].estimate)

        if dimension in affected_dimensions:
            drop_dimension = True
            for system in systems:
                if (not result[dimension][system].is_from_spherical_enrichment()) and (result[dimension][system].estimate == max_value):
                    drop_dimension = False
                    break

            if drop_dimension: 
                affected_dimensions.remove(dimension)


def print_results_for_graph_system_set(
        prime_powers,
        systems,
        lowest_coordinate_count,
        highest_coordinate_count,
        lowest_dimension,
        highest_dimension,
        fout,
):    
    print('Dealing with systems: {0}'.format(systems), file=sys.stderr)

    result = defaultdict(lambda : defaultdict(BorsukLowerEstimate))
    affected_dimensions = set()

    for n in range(lowest_coordinate_count, highest_coordinate_count + 1):        
        update_progress_bar(n - lowest_coordinate_count + 1, highest_coordinate_count - lowest_coordinate_count + 1)

        if SYSTEM_G1_THROW_OUT in systems:
            dimension = calculate_dimension_simple_dual_transition(n)
            result[dimension][SYSTEM_G1_THROW_OUT] = find_f_best_lower_estimate_fixed(n, prime_powers)
            affected_dimensions = try_update_affected_dimensions(result, dimension, SYSTEM_G1_THROW_OUT, affected_dimensions)

        if SYSTEM_G2_THROW_OUT in systems:
            dimension = calculate_dimension_simple_dual_transition(n)
            result[dimension][SYSTEM_G2_THROW_OUT] = find_f_best_lower_estimate_free(n, prime_powers)        
            affected_dimensions = try_update_affected_dimensions(result, dimension, SYSTEM_G2_THROW_OUT, affected_dimensions)
            
        remainder_modulo_four_group = get_remainder_modulo_four_group(n)
        if remainder_modulo_four_group == 0:
            if SYSTEM_MOD_0_FREE in systems:
                dimension = calculate_dimension_two_digits_single_ban_custom_remainder(n)
                result[dimension][SYSTEM_MOD_0_FREE] = find_f_two_digits_single_ban_custom_remainder(n, prime_powers)
                affected_dimensions = try_update_affected_dimensions(result, dimension, SYSTEM_MOD_0_FREE, affected_dimensions)
        elif remainder_modulo_four_group == 3:
            if SYSTEM_MOD_MINUS_1_FREE in systems:
                dimension = calculate_dimension_two_digits_single_ban_custom_remainder(n)
                result[dimension][SYSTEM_MOD_MINUS_1_FREE] = find_f_two_digits_single_ban_custom_remainder(n, prime_powers)
                affected_dimensions = try_update_affected_dimensions(result, dimension, SYSTEM_MOD_MINUS_1_FREE, affected_dimensions)
        elif remainder_modulo_four_group == 1:
            if SYSTEM_MOD_PLUS_1_FREE in systems:
                dimension = calculate_dimension_two_digits_single_ban_custom_remainder(n)
                result[dimension][SYSTEM_MOD_PLUS_1_FREE] = find_f_two_digits_single_ban_custom_remainder(n, prime_powers)
                affected_dimensions = try_update_affected_dimensions(result, dimension, SYSTEM_MOD_PLUS_1_FREE, affected_dimensions)
        else:
            assert remainder_modulo_four_group == 2

        for (system, dimension, estimate) in find_f_two_restrictions_common_remainder(n, prime_powers):
            if system in systems:
                result[dimension][system].update_if_better(estimate)
                affected_dimensions = try_update_affected_dimensions(result, dimension, system, affected_dimensions)

    perform_spherical_enrichment(1, highest_dimension, systems, result, affected_dimensions)
    make_bold(lowest_dimension, highest_dimension, systems, result)

    print_table(
        lowest_dimension, 
        highest_dimension, 
        systems, 
        result, 
        affected_dimensions,
        fout
    )


def main():
    parser = ArgumentParser()
    parser.add_argument('-p', '--prime-powers', help='File with consecutive prime powers, one prime power per line', default=DEFAULT_PRIME_POWERS_FILENAME)
    parser.add_argument('-lo', '--lowest', help='Lowest coordinate count', type=int, default=DEFAULT_LOWEST_COORDINATE_COUNT)
    parser.add_argument('-hi', '--highest', help='Highest coordinate count', type=int, default=DEFAULT_HIGHEST_COORDINATE_COUNT)
    parser.add_argument('-F', '--lowest-dimension', default=DEFAULT_LOWEST_DIMENSION, type=int)
    parser.add_argument('-T', '--highest-dimension', default=DEFAULT_HIGHEST_DIMENSION, type=int)
    parser.add_argument('-o', '--output', help='Output in LaTeX format', default='output.tex')
    args = parser.parse_args()

    prime_powers = load_prime_powers(args.prime_powers)

    with open(args.output, 'w') as fout:
        print(TEX_BEGIN, file=fout)
        print_results_for_graph_system_set(
            prime_powers,
            [SYSTEM_G1_THROW_OUT, SYSTEM_G2_THROW_OUT],
            args.lowest,
            args.highest,
            args.lowest_dimension,
            args.highest_dimension,
            fout,
        )
        print_results_for_graph_system_set(
            prime_powers,
            [SYSTEM_MOD_0_FREE, SYSTEM_MOD_MINUS_1_FREE, SYSTEM_MOD_PLUS_1_FREE],
            args.lowest,
            args.highest,
            args.lowest_dimension,
            args.highest_dimension,
            fout,
        )
        print_results_for_graph_system_set(
            prime_powers,
            [SYSTEM_TWO_RESTRICTIONS, SYSTEM_TWO_RESTRICTIONS_MINUS_1, SYSTEM_TWO_RESTRICTIONS_MINUS_2, SYSTEM_TWO_RESTRICTIONS_MINUS_3],
            args.lowest,
            args.highest,
            args.lowest_dimension,
            args.highest_dimension,
            fout,
        )
        print_results_for_graph_system_set(
            prime_powers,
            [SYSTEM_G2_THROW_OUT, SYSTEM_MOD_0_FREE, SYSTEM_MOD_PLUS_1_FREE, SYSTEM_TWO_RESTRICTIONS],
            args.lowest,
            args.highest,
            args.lowest_dimension,
            args.highest_dimension,
            fout,
        )
        print_results_for_graph_system_set(
            prime_powers,
            [SYSTEM_G2_THROW_OUT, SYSTEM_TWO_RESTRICTIONS],
            args.lowest,
            args.highest,
            args.lowest_dimension,
            args.highest_dimension,
            fout,
        )
        print(TEX_END, file=fout)


if __name__ == '__main__':
    main()
