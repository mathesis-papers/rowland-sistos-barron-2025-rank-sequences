import json
from math import comb


def multiply(left, right):
    product = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            product[left_degree + right_degree] += left_value * right_value
    return product


def multinomial_coefficients(parameter, exponent):
    coefficients = [0] * (2 * exponent + 1)
    for quadratic_count in range(exponent + 1):
        for linear_count in range(exponent - quadratic_count + 1):
            degree = 2 * quadratic_count + linear_count
            coefficients[degree] += (
                comb(exponent, quadratic_count)
                * comb(exponent - quadratic_count, linear_count)
                * (-parameter) ** quadratic_count
            )
    return coefficients


def sequence_value(parameter, base, index):
    return parameter * base ** (3 * index) - base ** (2 * index) - base**index


def determinant_bareiss(matrix):
    values = [row[:] for row in matrix]
    dimension = len(values)
    previous = 1
    sign = 1
    for column in range(dimension - 1):
        pivot_row = next(
            (row for row in range(column, dimension) if values[row][column]), None
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            values[column], values[pivot_row] = values[pivot_row], values[column]
            sign = -sign
        pivot = values[column][column]
        for row in range(column + 1, dimension):
            for next_column in range(column + 1, dimension):
                numerator = (
                    values[row][next_column] * pivot
                    - values[row][column] * values[column][next_column]
                )
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                values[row][next_column] = quotient
            values[row][column] = 0
        previous = pivot
    return sign * values[-1][-1]


def hankel_check(parameter, exponent, base):
    coefficients = multinomial_coefficients(parameter, exponent)
    active_degrees = [degree for degree, value in enumerate(coefficients) if value]
    dimension = len(active_degrees)
    values = [
        sequence_value(parameter, base, index) ** exponent
        for index in range(3 * dimension + 1)
    ]
    matrix = [
        [values[row + column] for column in range(dimension)]
        for row in range(dimension)
    ]
    determinant = determinant_bareiss(matrix)
    weights = [(-1) ** exponent * coefficients[degree] for degree in active_degrees]
    roots = [base ** (exponent + degree) for degree in active_degrees]
    expected_determinant = 1
    for weight in weights:
        expected_determinant *= weight
    for left_index, left_root in enumerate(roots):
        for right_root in roots[left_index + 1 :]:
            expected_determinant *= (right_root - left_root) ** 2
    assert determinant == expected_determinant != 0
    annihilator = [1]
    for root in roots:
        annihilator = multiply(annihilator, [-root, 1])
    for start in range(len(values) - dimension):
        assert sum(
            coefficient * values[start + offset]
            for offset, coefficient in enumerate(annihilator)
        ) == 0
    return {
        'parameter': parameter,
        'exponent': exponent,
        'base': base,
        'rank': dimension,
        'hankel_determinant_nonzero': True,
        'vandermonde_determinant_identity': True,
        'annihilating_recurrence_checks': len(values) - dimension,
    }


def general_products(initial_rank, spacing, exponent):
    dimension = initial_rank - 2
    zero = (0,) * dimension
    signatures = [zero, (1,) + zero[1:], (spacing,) + zero[1:]]
    for coordinate in range(1, dimension):
        signatures.append(tuple(int(index == coordinate) for index in range(dimension)))
    products = {zero: 1}
    for _ in range(exponent):
        expanded = {}
        for signature, coefficient in products.items():
            for root in signatures:
                result = tuple(left + right for left, right in zip(signature, root))
                expanded[result] = expanded.get(result, 0) + coefficient
        products = expanded
    return products


def general_family_checks():
    formula_checks = 0
    coefficient_sum_checks = 0
    evaluation_checks = 0
    first_collision_checks = 0
    examples = []
    prime_bases = [2, 3, 5, 7]
    for initial_rank in range(3, 7):
        for spacing in range(2, 13):
            for exponent in range(1, 13):
                products = general_products(initial_rank, spacing, exponent)
                expected = comb(exponent + initial_rank - 1, initial_rank - 1)
                if exponent >= spacing:
                    expected -= comb(exponent - spacing + initial_rank - 1, initial_rank - 1)
                assert len(products) == expected
                assert all(coefficient > 0 for coefficient in products.values())
                formula_checks += 1
                assert sum(products.values()) == initial_rank**exponent
                coefficient_sum_checks += 1
                if initial_rank <= 5 and spacing <= 5 and exponent <= 4:
                    for index in range(4):
                        original = 1 + 2**index + 2**(spacing * index)
                        original += sum(prime**index for prime in prime_bases[1:initial_rank-2])
                        expanded = 0
                        for signature, coefficient in products.items():
                            root = 1
                            for prime, valuation in zip(prime_bases, signature):
                                root *= prime**valuation
                            expanded += coefficient * root**index
                        assert expanded == original**exponent
                        evaluation_checks += 1
                if exponent == spacing:
                    assert expected == comb(exponent + initial_rank - 1, initial_rank - 1) - 1
                    assert len(general_products(initial_rank, spacing + 1, exponent)) == expected + 1
                    first_collision_checks += 1
                if initial_rank == 3 and spacing in [2, 3, 4, 5] and exponent == 12:
                    assert expected == spacing * exponent + spacing * (3 - spacing) // 2
                    examples.append({'initial_rank': initial_rank, 'spacing': spacing,
                                     'exponent': exponent, 'power_rank': expected})
    return {'domain': {'initial_rank': [3, 6], 'spacing': [2, 12], 'exponent': [1, 12]},
            'positive_product_support_formula_checks': formula_checks,
            'multinomial_coefficient_sum_checks': coefficient_sum_checks,
            'integer_sequence_expansion_checks': evaluation_checks,
            'first_collision_checks': first_collision_checks, 'examples': examples}


def main():
    coefficient_cases = 0
    coefficient_equalities = 0
    evaluation_equalities = 0
    cancellation_parameters = []
    prefixes = {}
    examples = []
    for parameter in range(3, 31):
        coefficients = [1]
        ranks = []
        for exponent in range(1, 66):
            coefficients = multiply(coefficients, [1, 1, -parameter])
            independent = multinomial_coefficients(parameter, exponent)
            assert coefficients == independent
            coefficient_cases += 1
            coefficient_equalities += len(coefficients)
            assert coefficients[0] == 1
            assert coefficients[1] == exponent
            assert coefficients[2] == comb(exponent, 2) - parameter * exponent
            rank = sum(value != 0 for value in coefficients)
            ranks.append(rank)
            if exponent == 2 * parameter + 1:
                assert coefficients[2] == 0
                assert rank <= 2 * exponent
                cancellation_parameters.append(parameter)
            if exponent <= 10:
                for base in [2, 3, 5]:
                    for index in range(5):
                        expanded = (-1) ** exponent * sum(
                            coefficient * base ** ((exponent + degree) * index)
                            for degree, coefficient in enumerate(coefficients)
                        )
                        assert expanded == sequence_value(parameter, base, index) ** exponent
                        evaluation_equalities += 1
            if parameter in [3, 11] and exponent == 7:
                examples.append({
                    'parameter': parameter,
                    'exponent': exponent,
                    'rank': rank,
                    'zero_degrees': [
                        degree for degree, value in enumerate(coefficients) if value == 0
                    ],
                })
        prefixes[str(parameter)] = ranks[:10]
    assert len(cancellation_parameters) == 28
    positivity_cases = 0
    for exponent in range(66):
        coefficients = multinomial_coefficients(-1, exponent)
        assert all(value > 0 for value in coefficients)
        positivity_cases += len(coefficients)
    hankel = [
        hankel_check(parameter, exponent, base)
        for parameter, exponent, base in [
            (3, 1, 2), (3, 2, 2), (3, 7, 2), (11, 7, 2), (3, 7, 3)
        ]
    ]
    result = {
        'schema_version': 1,
        'all_checks_passed': True,
        'coefficient_domain': {'parameter': [3, 30], 'exponent': [1, 65]},
        'coefficient_methods': ['iterated integer convolution', 'multinomial expansion'],
        'coefficient_polynomials_compared': coefficient_cases,
        'coefficient_equalities': coefficient_equalities,
        'exponential_expansion_equalities': evaluation_equalities,
        'exponential_expansion_domain': {
            'parameter': [3, 30], 'exponent': [1, 10], 'bases': [2, 3, 5], 'index': [0, 4]
        },
        'guaranteed_cancellation_parameters': cancellation_parameters,
        'positive_trinomial_coefficients_checked': positivity_cases,
        'examples': examples,
        'rank_prefixes_exponents_1_through_10': prefixes,
        'hankel_checks': hankel,
        'general_families': general_family_checks(),
        'scope': 'Finite coefficient identities, exponential expansions, exact Hankel determinants, and positive-product support counts; the all-parameter conclusions have written proofs in note.tex.'
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
