#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 08:18:33 2023

@author: greglong

# Inspired by the work of Luca Pasqualini (C) 2019: https://github.com/InsaneMonster/NistRng
# Inspired by the work of David Johnston (C) 2017:  https://github.com/dj-on-github/sp800_22_tests


"""

import numpy
import scipy
import math
import time
import copy
import sys
import random


def monobit(bits: numpy.ndarray):
    """
    Monobit test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of the test is the proportion of zeroes and ones for the entire sequence. The purpose of this test is to determine
    whether the number of ones and zeros in a sequence are approximately the same as would be expected for a truly random sequence.
    The test assesses the closeness of the fraction of ones to 1/2, that is, the number of ones and zeroes in a sequence
    should be about the same.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 100 bits (i.e., n ≥ 100). 
    
    It is recommended this is the first NIST randomness test to be ran and only proceed with other tests if this is passes
    
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()

    # Compute ones and zeroes
    ones: int = numpy.count_nonzero(bits)
    zeroes: int = bits.size - ones
    # Compute difference
    difference: int = abs(ones - zeroes)
    # Compute score
    score: float = math.erfc(float(difference) / (math.sqrt(float(bits.size)) * math.sqrt(2.0)))
    elapse_time = (time.time()-start_time)*1000
    # Return result
    if score >= significance_value:
        return ('Monobit', True, numpy.array(score),elapse_time)
    return ('Monobit', False, numpy.array(score),elapse_time) 

def FrequencyWithinBlock(bits: numpy.ndarray):
    """
    Frequency within block test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of the test is the proportion of ones within M-bit blocks. The purpose of this test is to determine whether the frequency of
    ones in an M-bit block is approximately M/2, as would be expected under an assumption of randomness.
    For block size M=1, this test degenerates to the Frequency (Monobit) test.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 100 bits (i.e., n ≥ 100).
    """    
    # Defining Test's Signficiant Value
    significance_value = 0.01
    
    # Define specific test attributes
    sequence_size_min: int = 100
    default_block_size: int = 20
    blocks_number_max: int = 100
    # Define cache attributes
    last_bits_size: int = -1
    block_size: int = -1
    blocks_number: int = -1
    
    # Defining Start Time
    start_time = time.time()

    if last_bits_size == -1 or last_bits_size != bits.size:
        # Get the number of blocks (N) with the default minimum block size (M)
        block_size: int = default_block_size
        blocks_number: int = int(bits.size // block_size)
        # Get the block size (M) if the number of blocks (N) exceed the allowed max
        if blocks_number >= blocks_number_max:
            blocks_number = blocks_number_max - 1
            block_size = int(bits.size // blocks_number)
        # Save in the cache
        last_bits_size = bits.size
        block_size = block_size
        blocks_number = blocks_number
    else:
        block_size: int = block_size
        blocks_number: int = blocks_number
    # Initialize a list of fractions
    block_fractions: numpy.ndarray = numpy.zeros(blocks_number, dtype=float)
    for i in range(blocks_number):
        # Get the bits in the current block
        block: numpy.ndarray = bits[i * block_size:((i + 1) * block_size)]
        # Compute ones and save the fraction in the array
        block_fractions[i] = numpy.count_nonzero(block) / block_size
    # Compute Chi-square
    chi_square: float = numpy.sum(4.0 * block_size * ((block_fractions[:] - 0.5) ** 2))
    # Compute score (P-value) applying the lower incomplete gamma function
    score: float = scipy.special.gammaincc((blocks_number / 2.0), chi_square / 2.0)
    elapse_time = (time.time()-start_time)*1000
    
    # Return result
    if score >= significance_value:
        return ('Frequency Block Test', True, numpy.array(score),elapse_time)
    return ('Frequency Block Test', False, numpy.array(score),elapse_time)

def Runs(bits: numpy.ndarray):
    """
    Runs test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the total number of runs in the sequence, where a run is an uninterrupted sequence of identical bits.
    A run of length k consists of exactly k identical bits and is bounded before and after with a bit of the opposite value.
    The purpose of the runs test is to determine whether the number of runs of ones and zeros of various lengths is as expected
    for a random sequence. In particular, this test determines whether the oscillation between such zeros and ones is too fast or too slow.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 100 bits (i.e., n ≥ 100). 
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    proportion: float = numpy.count_nonzero(bits) / bits.size
    # Count the observed runs (list of adjacent equal bits)
    observed_runs: float = 1.0
    for i in range(bits.size - 1):
        if bits[i] != bits[i + 1]:
            observed_runs += 1.0
    # Compute score (P-value)
    score: float = math.erfc(abs(observed_runs - (2.0 * bits.size * proportion * (1.0 - proportion))) / (2.0 * math.sqrt(2.0 * bits.size) * proportion * (1 - proportion)))
    elapse_time = (time.time()-start_time)*1000
    
    # Return result
    if score >= significance_value:
        return ('Runs', True, numpy.array(score), elapse_time)
    return ('Runs', False, numpy.array(score), elapse_time)

def _probabilities(size_of_block: int, index: int) -> float:
    """
    Returns a probability at the given index in the array or probabilities defined for the block of the given size.

    :param size_of_block: can be 8, 128, 512, 1000 and in any other case will fallback on 10000
    :param index: the index of the probability
    :return: the probability at the given index
    """
    if size_of_block == 8:
        return [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124][index]
    elif size_of_block == 128:
        return [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124][index]
    elif size_of_block == 512:
        return [0.1170, 0.2460, 0.2523, 0.1755, 0.1027, 0.1124][index]
    elif size_of_block == 1000:
        return [0.1307, 0.2437, 0.2452, 0.1714, 0.1002, 0.1088][index]
    else:
        return [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727][index]

def OnesInABlock(bits: numpy.ndarray):
    """
    Longest run ones in a block test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of the test is the longest run of ones within M-bit blocks. The purpose of this test is to determine whether
    the length of the longest run of ones within the tested sequence is consistent with the length of the longest run of
    ones that would be expected in a random sequence. Note that an irregularity in the expected length of the longest run
    of ones implies that there is also an irregularity in the expected length of the longest run of zeroes.
    Therefore, only a test for ones is necessary.

    The significance value of the test is 0.01.
    
    This test has predefined blocks when N is of a certain size.  Recommended N > 6272 for best results
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Define specific test attributes
    sequence_size_min: int = 128
    # Define cache attributes
    last_bits_size: int = -1
    block_size: int = -1
    blocks_number: int = -1
    k: int = -1
    # Defining Start Time
    start_time = time.time()
    
    # Reload values is cache is empty or no longer up-to-date
    # Otherwise, use cache
    if last_bits_size == -1 or last_bits_size != bits.size:
        # Set the block size depending on the inumpyut sequence length
        block_size: int = 10000
        if bits.size < 6272:
            block_size: int = 8
        elif bits.size < 750000:
            block_size: int = 128
        # Set the block number and K depending on the block size
        k: int = 6
        blocks_number: int = 75
        if block_size == 8:
            k: int = 3
            blocks_number: int = 16
        elif block_size == 128:
            k: int = 5
            blocks_number: int = 49
        # Save in the cache
        last_bits_size = bits.size
        block_size = block_size
        blocks_number = blocks_number
        k = k
    else:
        block_size: int = block_size
        blocks_number: int = blocks_number
        k: int = k
    # Define the array of frequencies
    frequencies: numpy.ndarray = numpy.zeros(7, dtype=int)
    # Find longest run length in each block
    for i in range(blocks_number):
        block: numpy.ndarray = bits[i * block_size:((i + 1) * block_size)]
        run_length: int = 0
        longest_run_length: int = 0
        # Count the length of each adjacent bits group (runs) in the current block and update the max length of them
        for j in range(block_size):
            if block[j] == 1:
                run_length += 1
                if run_length > longest_run_length:
                    longest_run_length = run_length
            else:
                run_length = 0
        # Update the list of frequencies
        if block_size == 8:
            frequencies[min(3, max(0, longest_run_length - 1))] += 1
        elif block_size == 128:
            frequencies[min(5, max(0, longest_run_length - 4))] += 1
        else:
            frequencies[min(6, max(0, longest_run_length - 10))] += 1
    # Compute Chi-square
    chi_square: float = 0.0
    for i in range(k + 1):
        chi_square += ((frequencies[i] - blocks_number * _probabilities(block_size, i)) ** 2) / (blocks_number * _probabilities(block_size, i))
    # Compute score (P-value)
    score: float = scipy.special.gammaincc(k / 2.0, chi_square / 2.0)
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Ones in a Run', True, numpy.array(score), elapse_time)
    return ('Ones in a Run', False, numpy.array(score), elapse_time)


def perform_row_operations(matrix, i, forward_elimination):
    """
    Performs elementary row operations.
    """
    if forward_elimination:
        j = i + 1
        while j < matrix.shape[0]:
            if matrix[j][i] == 1:
                matrix[j, :] = (matrix[j, :] + matrix[i, :]) % 2
            j += 1
    else:
        j = i - 1
        while j >= 0:
            if matrix[j][i] == 1:
                matrix[j, :] = (matrix[j, :] + matrix[i, :]) % 2
            j -= 1

def find_unit_element_swap(matrix, i, forward_elimination):
    """
    Searches through the rows below/above the given index to see which rows contain 1.
    """
    row_swap_operation = 0
    if forward_elimination:
        index = i + 1
        while index < matrix.shape[0] and matrix[index][i] == 0:
            index += 1
        if index < matrix.shape[0]:
            row_swap_operation = swap_rows(matrix, i, index)
    else:
        index = i - 1
        while index >= 0 and matrix[index][i] == 0:
            index -= 1
        if index >= 0:
            row_swap_operation = swap_rows(matrix, i, index)
    return row_swap_operation

def swap_rows(matrix, source_row_index, target_row_index):
    """
    Swaps two rows in a matrix.
    """
    temp_matrix = copy.copy(matrix[source_row_index, :])
    matrix[source_row_index, :] = matrix[target_row_index, :]
    matrix[target_row_index, :] = temp_matrix
    return 1

def compute_rank(matrix, base_rank):
    """
    Computes the rank of the transformed matrix.
    """
    rank = base_rank
    i = 0
    while i < matrix.shape[0]:
        all_zeros = all(matrix[i] == 0)
        if all_zeros:
            rank -= 1
        i += 1
    return rank

def product(number_of_rows, number_of_cols):
    indexes = numpy.arange(number_of_rows)
    product_value = float(numpy.prod(((1.0 - (2.0 ** (indexes[:] - number_of_cols))) * (1.0 - (2.0 ** (indexes[:] - number_of_rows)))) / (1 - (2.0 ** (indexes[:] - number_of_rows)))))
    return product_value

def binary_matrix_rank(bits: numpy.ndarray):
    
    """
    The focus of the test is the rank of disjoint sub-matrices of the entire sequence. The purpose of this test is
    to check for linear dependence among fixed length substrings of the original sequence. Note that this test
    also appears in the DIEHARD battery of tests [7]. 
    
    Significance value of this test is set at 0.01
    
    For rows_number = cols_number = 32, each sequence to be tested should consist of a minimum of 38,912 bits. 
    """
    
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    rows_number = 32
    cols_number = 32
    block_size_min = 38
    full_rank_probability = product(rows_number, cols_number) * (2.0 ** ((rows_number * (cols_number + rows_number - rows_number)) - (rows_number * cols_number)))
    minus_rank_probability = product(rows_number - 1, cols_number) * (2.0 ** ((rows_number * (cols_number + rows_number - rows_number)) - (rows_number * cols_number)))
    remained_rank_probability = 1.0 - (full_rank_probability + minus_rank_probability)

    last_bits_size = -1
    blocks_number = -1

    if last_bits_size == -1 or last_bits_size != bits.size:
        blocks_number = int(math.floor(bits.size / (rows_number * cols_number)))
        last_bits_size = bits.size
    else:
        blocks_number = blocks_number

    full_rank_matrices = 0
    minus_rank_matrices = 0
    remainder = 0

    for i in range(blocks_number):
        block = bits[i * (rows_number * cols_number):(i + 1) * (rows_number * cols_number)].reshape((rows_number, cols_number))
        matrix = block
        base_rank = min(rows_number, cols_number)
        i = 0
        while i < base_rank - 1:
            if matrix[i][i] == 1:
                perform_row_operations(matrix, i, True)
            else:
                found = find_unit_element_swap(matrix, i, True)
                if found == 1:
                    perform_row_operations(matrix, i, True)
            i += 1
        i = base_rank - 1
        while i > 0:
            if matrix[i][i] == 1:
                perform_row_operations(matrix, i, False)
            else:
                if find_unit_element_swap(matrix, i, False) == 1:
                    perform_row_operations(matrix, i, False)
            i -= 1
        rank = compute_rank(matrix, base_rank)

        if rank == rows_number:
            full_rank_matrices += 1
        elif rank == rows_number - 1:
            minus_rank_matrices += 1
        else:
            remainder += 1

    chi_square = (((full_rank_matrices - (full_rank_probability * blocks_number)) ** 2) / (full_rank_probability * blocks_number)) + \
                 (((minus_rank_matrices - (minus_rank_probability * blocks_number)) ** 2) / (minus_rank_probability * blocks_number)) + \
                 (((remainder - (remained_rank_probability * blocks_number)) ** 2) / (remained_rank_probability * blocks_number))

    score = math.e ** (-chi_square / 2.0)
    elapse_time = (time.time()-start_time)*1000

    if score >= significance_value:
        return ('Binary Matrix Rank', True, numpy.array(score), elapse_time)
    return ('Binary Matrix Rank', False, numpy.array(score), elapse_time)

def dft(bits: numpy.ndarray):
    """
    Discrete Fourier transform (spectral) test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the peak heights in the Discrete Fourier Transform of the sequence.
    The purpose of this test is to detect periodic features (i.e., repetitive patterns that are near each other) in the
    tested sequence that would indicate a deviation from the assumption of randomness.
    The intention is to detect whether the number of peaks exceeding the 95% threshold is significantly different than 5%.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 1000 bits (i.e., n ≥ 1000). 
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    # Make sure the sequence is even in length
    bits_copy: numpy.ndarray = bits.copy()
    if (bits_copy.size % 2) == 1:
        bits_copy = bits_copy[:-1]
    # Convert all the zeros in the array to -1
    bits_copy[bits_copy == 0] = -1
    # Compute DFT
    discrete_fourier_transform = numpy.fft.fft(bits_copy)
    # Compute magnitudes of first half of sequence depending on the system type
    if sys.version_info > (3, 0):
        magnitudes = abs(discrete_fourier_transform)[:bits_copy.size // 2]
    else:
        magnitudes = abs(discrete_fourier_transform)[:bits_copy.size / 2]
    # Compute upper threshold
    threshold: float = math.sqrt(math.log(1.0 / 0.05) * bits_copy.size)
    # Compute the expected number of peaks (N0)
    expected_peaks: float = 0.95 * bits_copy.size / 2.0
    # Count the peaks above the upper threshold (N1)
    counted_peaks: float = float(len(magnitudes[magnitudes < threshold]))
    # Compute the score (P-value) using the normalized difference
    normalized_difference: float = (counted_peaks - expected_peaks) / math.sqrt((bits_copy.size * 0.95 * 0.05) / 4)
    score: float = math.erfc(abs(normalized_difference) / math.sqrt(2))
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Discrete Fourier Transform', True, numpy.array(score), elapse_time)
    return ('Discrete Fourier Transform', False, numpy.array(score), elapse_time)


def choose_template(templates):
    return numpy.array(random.choice(random.choice(templates)))

def split_bits_into_blocks(bits, blocks_number):
    substring_bits_length = int(bits.size // blocks_number)
    return bits.size, substring_bits_length

def count_matches(block, b_template):
    position, count = 0, 0
    while position < (len(block) - len(b_template)):
        if (block[position:position + len(b_template)] == b_template).all():
            position += len(b_template)
            count += 1
        else:
            position += 1
    return count

def compute_mu_sigma(substring_bits_length, b_template_size):
    mu = float(substring_bits_length - b_template_size + 1) / float(2 ** b_template_size)
    sigma = substring_bits_length * ((1.0 / float(2 ** b_template_size)) - (float((2 * b_template_size) - 1) / float(2 ** (2 * b_template_size))))
    return mu, sigma

def compute_chi_square(matches, mu, sigma):
    return float(numpy.sum(((matches - mu) ** 2) / (sigma ** 2)))

def compute_score(chi_square, blocks_number):
    return scipy.special.gammaincc(blocks_number / 2.0, chi_square / 2.0)

def non_overlapping_match(bits: numpy.ndarray):
    """
    Non overlapping template matching test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the number of occurrences of  pre-specified target strings.
    The purpose of this test is to detect generators that produce too many occurrences of a given non-periodic (aperiodic)
    pattern. For this test an m-bit window is used to search for a specific m-bit pattern. If the pattern is not found,
    the window slides one bit position. If the pattern is found, the window is reset to the bit after the found pattern,
    and the search resumes.

    The significance value of the test is 0.01.
    
    Note: at least 1,028,016 bits are required
    """
    
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    templates = [
        [[0, 1], [1, 0]],
        [[0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0]],
        # ... other templates ...
    ]
    blocks_number = 8
    significance_value = 0.01

    b_template = choose_template(templates)
    last_bits_size, substring_bits_length = split_bits_into_blocks(bits, blocks_number)

    matches = numpy.zeros(blocks_number, dtype=int)

    for i in range(blocks_number):
        block = bits[i * substring_bits_length:(i + 1) * substring_bits_length]
        matches[i] = count_matches(block, b_template)

    mu, sigma = compute_mu_sigma(substring_bits_length, len(b_template))
    chi_square = compute_chi_square(matches, mu, sigma)

    if chi_square != 0:
        score = compute_score(chi_square, blocks_number)
    
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Non-Overlaping Match', True, numpy.array(score), elapse_time)
    return ('Non-Overlaping Match', False, numpy.array(score), elapse_time)

def _log_gamma(x: []) -> []:
    return numpy.log(scipy.special.gamma(x))

def _get_probabilities(freedom_degree_values: [], eta_value: float) -> []:
    probabilities = []
    for freedom_degree_value in freedom_degree_values:
        if freedom_degree_value == 0:
            probability = numpy.exp(-eta_value)
        else:
            indexes = numpy.arange(1, freedom_degree_value + 1)
            probability = float(numpy.sum(numpy.exp(-eta_value - freedom_degree_value * numpy.log(2) + indexes[:] * numpy.log(eta_value) - _log_gamma(indexes[:] + 1) + _log_gamma(freedom_degree_value) - _log_gamma(indexes[:]) - _log_gamma(freedom_degree_value - indexes[:] + 1))))
        probabilities.append(probability)
    return probabilities

def overlapping_matching(bits: numpy.ndarray):
    """
    Overlapping template matching test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of the Overlapping Template Matching test is the number of occurrences of pre-specified target strings.
    Both this test and the Non-overlapping Template Matching test use an m-bit window to search for a specific m-bit pattern.
    As with the other test, if the pattern is not found, the window slides one bit position.
    The difference between this test and the other is that when the pattern is found, the window slides only one bit before resuming the search.

    The significance value of the test is 0.01.
    
    The values of K, M and N have been chosen such that each sequence to be tested consists of a minimum
    of 10^6 bits (i.e., n ≥ 1,000,000)
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    template_bits_length = 10
    blocks_number = 968
    freedom_degrees = 5
    substring_bits_length = 1062

    # Build the template B as a fixed-sized sequence of ones
    b_template = numpy.ones(template_bits_length, dtype=int)
    # Count the distribution of matches of the template across blocks
    matches_distributions = numpy.zeros(freedom_degrees + 1, dtype=int)
    
    for i in range(blocks_number):
        # Define the block at the current index
        block = bits[i * substring_bits_length:(i + 1) * substring_bits_length]
        # Define counting variable
        count = 0
        # Count the matches in the block with respect to the given template
        for position in range(substring_bits_length - template_bits_length):
            if numpy.all(block[position:position + template_bits_length] == b_template):
                count += 1
        matches_distributions[min(count, freedom_degrees)] += 1

    # Define eta and default probabilities (from STS) of size freedom degrees + 1
    eta = (substring_bits_length - template_bits_length + 1.0) / (2.0 ** template_bits_length) / 2.0
    probabilities = numpy.array([0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865])
    
    # Compute probabilities up to degrees of freedom and change the last based on the sum of all of them
    probabilities[:freedom_degrees] = _get_probabilities(numpy.arange(freedom_degrees)[:], eta)
    probabilities[-1] = 1.0 - numpy.sum(probabilities)

    # Compute Chi-square
    chi_square = float(numpy.sum(((matches_distributions[:] - (blocks_number * probabilities[:])) ** 2) / (blocks_number * probabilities[:])))
    

    if chi_square != 0:
        # Compute the score (P-value)
        score: float = scipy.special.gammaincc(5.0 / 2.0, chi_square / 2.0)
        # Return result
        elapse_time = (time.time()-start_time)*1000

        if score >= significance_value:
            return ('Overlaping Match', True, numpy.array(score), elapse_time)
    return ('Overlaping Match', False, numpy.array(0), elapse_time)

def pattern_to_int(bit_pattern: numpy.ndarray):
    result = 0
    for bit in bit_pattern:
        result = (result << 1) + bit
    return result

def maurers_universal(bits: numpy.ndarray):
    """
    Maurers universal test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the number of bits between matching patterns (a measure that is related to the length of a compressed sequence).
    The purpose of the test is to detect whether or not the sequence can be significantly compressed without loss of information.
    A significantly compressible sequence is considered to be non-random.

    The significance value of the test is 0.01.
    
    The recommended size for this test is > 387,240 bits
    """
    
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    sequence_size_min = 387840
    default_pattern_size = 6
    thresholds = [904960, 2068480, 4654080, 10342400, 22753280, 49643520, 107560960, 231669760, 496435200, 1059061760]
    expected_value_table = [0, 0.73264948, 1.5374383, 2.40160681, 3.31122472, 4.25342659, 5.2177052, 6.1962507, 7.1836656, 8.1764248, 9.1723243, 10.170032, 11.168765, 12.168070, 13.167693, 14.167488, 15.167379]
    variance_table = [0, 0.690, 1.338, 1.901, 2.358, 2.705, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384, 3.401, 3.410, 3.416, 3.419, 3.421]

    if bits.size < sequence_size_min:
        return 0.0

    pattern_length = default_pattern_size
    for threshold in thresholds:
        if bits.size >= threshold:
            pattern_length += 1

    blocks_number = int(bits.size // pattern_length)
    q_blocks = 10 * (2 ** pattern_length)
    k_blocks = blocks_number - q_blocks

    table = numpy.zeros(2 ** pattern_length, dtype=int)
    computed_sum = 0.0

    for i in range(q_blocks):
        pattern = bits[i * pattern_length:(i + 1) * pattern_length]
        table[pattern_to_int(pattern)] = i + 1

    for i in range(q_blocks, blocks_number):
        pattern = bits[i * pattern_length:(i + 1) * pattern_length]
        difference = i + 1 - table[pattern_to_int(pattern)]
        table[pattern_to_int(pattern)] = i + 1
        computed_sum += math.log(difference, 2)

    fn = computed_sum / k_blocks
    magnitude = abs((fn - expected_value_table[pattern_length]) / ((math.sqrt(variance_table[pattern_length])) * math.sqrt(2)))
    score = math.erfc(magnitude)

    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Mauers Universal', True, numpy.array(score), elapse_time)
    return ('Mauers Universal', False, numpy.array(score), elapse_time)

def berlekamp_massey(sequence: numpy.ndarray):
    """
    Compute the linear complexity of a sequence of bits by the means of the Berlekamp Massey algorithm.

    :param sequence: the sequence of bits to compute the linear complexity for
    :return: the int value of the linear complexity
    """
    # Initialize b and c to all zeroes with first element one
    b = numpy.zeros(sequence.size, dtype=int)
    c = numpy.zeros(sequence.size, dtype=int)
    b[0] = 1
    c[0] = 1

    # Initialize the generator length
    generator_length = 0

    # Initialize variables
    m = -1
    n = 0
    while n < sequence.size:
        # Compute discrepancy
        discrepancy = sequence[n]
        for j in range(1, generator_length + 1):
            discrepancy = discrepancy ^ (c[j] & sequence[n - j])

        # If discrepancy is not zero, adjust polynomial
        if discrepancy != 0:
            t = c[:]
            for j in range(0, sequence.size - n + m):
                c[n - m + j] = c[n - m + j] ^ b[j]
            if generator_length <= n / 2:
                generator_length = n + 1 - generator_length
                m = n
                b = t
        n = n + 1

    # Return the length of generator
    return generator_length

def linear_complexity(bits: numpy.ndarray):
    """
    Linear complexity test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the length of a linear feedback shift register (LFSR). The purpose of this test is to determine whether
    or not the sequence is complex enough to be considered random. Random sequences are characterized by longer LFSRs.
    An LFSR that is too short implies non-randomness.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 1,000,000 bits (i.e., n ≥ 1,000,000).
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    # Define specific test attributes
    sequence_size_min = 1000000
    pattern_length = 512
    freedom_degrees = 6
    probabilities = numpy.array([0.010417, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833])
    # Compute mean
    mu = (pattern_length / 2.0) + (((-1) ** (pattern_length + 1)) + 9.0) / 36.0 - ((pattern_length / 3.0) + (2.0 / 9.0)) / (2 ** pattern_length)
    
    # Reload values if cache is empty or no longer up-to-date
    # Otherwise, use cache
    last_bits_size = -1
    blocks_number = -1
    if last_bits_size == -1 or last_bits_size != bits.size:
        # Set the block size
        blocks_number = int(bits.size // pattern_length)
        # Save in the cache
        last_bits_size = bits.size
        blocks_number = blocks_number
    else:
        blocks_number = blocks_number

    # Compute the linear complexity of the blocks
    blocks_linear_complexity = numpy.zeros(blocks_number, dtype=int)
    for i in range(blocks_number):
        blocks_linear_complexity[i] = berlekamp_massey(bits[(i * pattern_length):((i + 1) * pattern_length)])

    # Count the distribution over tickets
    tickets = ((-1.0) ** pattern_length) * (blocks_linear_complexity[:] - mu) + (2.0 / 9.0)

    # Compute frequencies depending on tickets
    frequencies = numpy.zeros(freedom_degrees + 1, dtype=int)
    for ticket in tickets:
        frequencies[min(freedom_degrees, int(max(-2.5, ticket) + 2.5))] += 1

    # Compute Chi-square using pre-defined probabilities
    chi_square = float(numpy.sum(((frequencies[:] - (blocks_number * probabilities[:])) ** 2.0) / (blocks_number * probabilities[:])))

    # Compute the score (P-value)
    score = scipy.special.gammaincc((freedom_degrees / 2.0), (chi_square / 2.0))

    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Linear Complexity', True, numpy.array(score), elapse_time)
    return ('Linear Complexity', False, numpy.array(score), elapse_time)

def count_pattern(pattern: numpy.ndarray, padded_sequence: numpy.ndarray, sequence_size: int) -> int:
    count = 0
    for i in range(sequence_size):
        match = all(pattern[j] == padded_sequence[i + j] for j in range(len(pattern)))
        if match:
            count += 1
    return count

def psi_sq_mv1(block_size: int, sequence_size: int, padded_sequence: numpy.ndarray) -> float:
    # Count the patterns
    counts = numpy.zeros(2 ** block_size, dtype=int)
    for i in range(2 ** block_size):
        pattern = (i >> numpy.arange(block_size, dtype=int)) & 1
        counts[i] = count_pattern(pattern, padded_sequence, sequence_size)

    # Compute Psi-Squared statistics and return it
    psi_sq_m = numpy.sum(counts[:] ** 2)
    psi_sq_m *= (2 ** block_size) / sequence_size
    psi_sq_m -= sequence_size
    return psi_sq_m

def serial(bits: numpy.ndarray):
    """
    Serial test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the frequency of all possible overlapping m-bit patterns across the entire sequence.
    The purpose of this test is to determine whether the number of occurrences of the 2mm-bit overlapping patterns is
    approximately the same as would be expected for a random sequence. Random sequences have uniformity; that is, every m-bit
    pattern has the same chance of appearing as every other m-bit pattern.
    Note that for m = 1, the Serial test is equivalent to the Monobit test.
    
    The significance value of the test is 0.01.
    
    Choose blocks_length and pattern_length such that blocks_length < [log2 n] -2. 
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    blocks_length_min = 4
    pattern_length = 4

    # Pad the sequence
    padded_bits = numpy.concatenate((bits, bits[0:pattern_length - 1]))

    # Compute Psi-Squared statistics
    psi_sq_m_0 = psi_sq_mv1(pattern_length, bits.size, padded_bits)
    psi_sq_m_1 = psi_sq_mv1(pattern_length - 1, bits.size, padded_bits)
    psi_sq_m_2 = psi_sq_mv1(pattern_length - 2, bits.size, padded_bits)
    delta_1 = psi_sq_m_0 - psi_sq_m_1
    delta_2 = psi_sq_m_0 - (2 * psi_sq_m_1) + psi_sq_m_2

    # Compute the scores (P-values)
    score_1 = scipy.special.gammaincc(2 ** (pattern_length - 2), delta_1 / 2.0)
    score_2 = scipy.special.gammaincc(2 ** (pattern_length - 3), delta_2 / 2.0)

    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score_1 >= significance_value and score_2 >= significance_value:
        return ('Serial', True, numpy.array([score_1]), elapse_time)
    return ('Serial', False, numpy.array([score_1]), elapse_time)

def _pattern_to_int(bit_pattern: numpy.ndarray) -> int:
    """
    Convert the given pattern of bits to an integer value.

    :param bit_pattern: the bit pattern to convert
    :return: the integer value identifying the pattern
    """
    result = 0
    for bit in bit_pattern:
        result = (result << 1) + bit
    return result

def approximate_entropy(bits: numpy.ndarray) -> tuple:
    """
    Approximate entropy test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    As with the Serial test, the focus of this test is the frequency of all possible overlapping m-bit patterns across the entire sequence.
    The purpose of the test is to compare the frequency of overlapping blocks of two consecutive/adjacent lengths (m and m+1) against the
    expected result for a random sequence.

    The significance value of the test is 0.01.
    
    Choose blocks_length such that blocks_length  < [log2 n]-5. 
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    
    blocks_length_min = 4
    pattern_length = 4

    # Reload values if cache is empty or no longer up-to-date
    # Otherwise, use cache
    last_bits_size = -1
    blocks_length = -1
    if last_bits_size == -1 or last_bits_size != bits.size:
        # Define the block length in a range bounded by 2 and 3
        blocks_length = min(2, max(3, int(math.floor(math.log(bits.size, 2))) - 6))
        # Save in the cache
        last_bits_size = bits.size
    else:
        blocks_length = blocks_length

    # Define Phi-m statistics list
    phi_m = []
    for iteration in range(blocks_length, blocks_length + 2):
        # Compute the padded sequence of bits
        padded_bits = numpy.concatenate((bits, bits[0:iteration - 1]))
        # Compute the frequency count
        counts = numpy.zeros(2 ** iteration, dtype=int)
        for i in range(2 ** iteration):
            count = sum(1 for j in range(bits.size) if _pattern_to_int(padded_bits[j:j + iteration]) == i)
            counts[i] = count
        # Compute C-i as the average of counts on the number of bits
        c_i = counts / float(bits.size)
        # Compute Phi-m based on C-i
        phi_m.append(numpy.sum(c_i[c_i > 0.0] * numpy.log((c_i[c_i > 0.0] / 10.0))))

    # Compute Chi-Square from the computed statistics
    chi_square = 2 * bits.size * (math.log(2) - (phi_m[0] - phi_m[1]))

    # Compute the score (P-value)
    score = scipy.special.gammaincc(2 ** (blocks_length - 1), (chi_square / 2.0))

    # Return result
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score >= significance_value:
        return ('Entropy Test', True, numpy.array(score), elapse_time)
    return ('Entropy', False, numpy.array(score), elapse_time)

def _compute_p_value(sequence_size: int, max_excursion: int) -> float:
        """
        Compute P-Value given the sequence size and the max excursion.

        :param sequence_size: the length of the sequence of bits
        :param max_excursion: the max excursion backward or forward
        :return: the computed float P-Value
        """
        # Execute first sum
        sum_a: float = 0.0
        start_k: int = int(math.floor((((float(-sequence_size) / max_excursion) + 1.0) / 4.0)))
        end_k: int = int(math.floor((((float(sequence_size) / max_excursion) - 1.0) / 4.0)))
        for k in range(start_k, end_k + 1):
            c: float = 0.5 * math.erfc(-(((4.0 * k) + 1.0) * max_excursion) / math.sqrt(sequence_size) * math.sqrt(0.5))
            d: float = 0.5 * math.erfc(-(((4.0 * k) - 1.0) * max_excursion) / math.sqrt(sequence_size) * math.sqrt(0.5))
            sum_a = sum_a + c - d
        # Execute second sum
        sum_b: float = 0.0
        start_k = int(math.floor((((float(-sequence_size) / max_excursion) - 3.0) / 4.0)))
        end_k = int(math.floor((((float(sequence_size) / max_excursion) - 1.0) / 4.0)))
        for k in range(start_k, end_k + 1):
            c: float = 0.5 * math.erfc(-(((4.0 * k) + 3.0) * max_excursion) / math.sqrt(sequence_size) * math.sqrt(0.5))
            d: float = 0.5 * math.erfc(-(((4.0 * k) + 1.0) * max_excursion) / math.sqrt(sequence_size) * math.sqrt(0.5))
            sum_b = sum_b + c - d
        # Return value
        return 1.0 - sum_a + sum_b

def cumulative_sum(bits: numpy.ndarray):
    """
    Cumulative sums test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the maximal excursion (from zero) of the random walk defined by the cumulative sum of adjusted (-1, +1) digits in the sequence.
    The purpose of the test is to determine whether the cumulative sum of the partial sequences occurring in the tested sequence is too large or too small
    relative to the expected behavior of that cumulative sum for random sequences.
    This cumulative sum may be considered as a random walk. For a random sequence, the excursions of the random walk should
    be near zero. For certain types of non-random sequences, the excursions of this random walk from zero will be large.
    
    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 100 bits (i.e., n ≥ 100). 
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    # Copy the bits to a new array
    bits_copy: numpy.ndarray = bits.copy()
    # Convert all the zeros in the array to -1
    bits_copy[bits_copy == 0] = -1
    # Compute the partial sum with forward (mode 0) and backward (mode 1) modes and record the largest excursion
    forward_sum: int = 0
    backward_sum: int = 0
    forward_max: int = 0
    backward_max: int = 0
    for i in range(bits_copy.size):
        forward_sum += bits_copy[i]
        backward_sum += bits_copy[bits_copy.size - 1 - i]
        forward_max = max(abs(forward_sum), forward_max)
        backward_max = max(abs(backward_sum), backward_max)
    # Compute the scores (P-Values)
    score_1: float = _compute_p_value(bits_copy.size, forward_max)
    score_2: float = _compute_p_value(bits_copy.size, backward_max)
    # Return result
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if score_1 >= significance_value and score_2 >= significance_value:
        return ('Cumulative Sum', True, numpy.array([score_1, score_2]), elapse_time)
    return ('Cumulative Sum', False, numpy.array([score_1, score_2]), elapse_time)

def random_excursion(bits: numpy.ndarray):
    """
    Random excursion test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the number of cycles having exactly K visits in a cumulative sum random walk.
    The cumulative sum random walk is derived from partial sums after the (0,1) sequence is transferred to the appropriate (-1, +1) sequence.
    A cycle of a random walk consists of a sequence of steps of unit length taken at random that begin at and return to the origin.
    The purpose of this test is to determine if the number of visits to a particular state within a cycle deviates from what one would expect
    for a random sequence. This test is actually a series of eight tests (and conclusions), one test and conclusion for each of the
    states: -4, -3, -2, -1 and +1, +2, +3, +4.
    
    The significance value of the test is 0.01.
    
    Recommended input size: 1,000,000 bits
    """
    _probabilities_xk = [numpy.array([0.5, 0.25, 0.125, 0.0625, 0.0312, 0.0312]),
                                  numpy.array([0.75, 0.0625, 0.0469, 0.0352, 0.0264, 0.0791]),
                                  numpy.array([0.8333, 0.0278, 0.0231, 0.0193, 0.0161, 0.0804]),
                                  numpy.array([0.875, 0.0156, 0.0137, 0.012, 0.0105, 0.0733]),
                                  numpy.array([0.9, 0.01, 0.009, 0.0081, 0.0073, 0.0656]),
                                  numpy.array([0.9167, 0.0069, 0.0064, 0.0058, 0.0053, 0.0588]),
                                  numpy.array([0.9286, 0.0051, 0.0047, 0.0044, 0.0041, 0.0531])]
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()

    #Copy the bits to a new array
    bits_copy: numpy.ndarray = bits.copy()
    # Convert all the zeros in the array to -1
    bits_copy[bits_copy == 0] = -1
    # Generate the padded cumulative sum of the array of -1, 1
    sum_prime: numpy.ndarray = numpy.concatenate((numpy.array([0]), numpy.cumsum(bits_copy), numpy.array([0]))).astype(int)
    # Compute the cycles iterating over each position of S' (sum_prime) and define the first cycle
    cycles: [] = []
    cycle: [] = [0]
    for index, _ in enumerate(sum_prime[1:]):
        # Once a zero crossing is found add all the non zero elements of S' to the cycle
        # Else wrap up the cycle and start a new cycle
        if sum_prime[index] != 0:
            cycle += [sum_prime[index]]
        else:
            cycle += [0]
            cycles.append(cycle)
            cycle: [] = [0]
    # Append the last cycle
    cycles.append(cycle)
    # Compute the size of the cycles list
    cycles_size: int = len(cycles)
    # Setup frequencies table (Vk(x))
    frequencies_table: dict = {
        -4: numpy.zeros(6, dtype=int),
        -3: numpy.zeros(6, dtype=int),
        -2: numpy.zeros(6, dtype=int),
        -1: numpy.zeros(6, dtype=int),
        1: numpy.zeros(6, dtype=int),
        2: numpy.zeros(6, dtype=int),
        3: numpy.zeros(6, dtype=int),
        4: numpy.zeros(6, dtype=int),
    }
    # Count occurrences
    for value in frequencies_table.keys():
        for k in range(frequencies_table[value].size):
            count: int = 0
            # Count how many cycles in which x occurs k times
            for cycle in cycles:
                # Count how many times the value used as key of the table occurs in the current cycle
                occurrences: int = numpy.count_nonzero(numpy.array(cycle) == value)
                # If the value occurs k times, increment the cycle count
                if 5 > k == occurrences:
                    count += 1
                elif occurrences >= 5:
                    count += 1
            frequencies_table[value][k] = count
    # Compute the scores (P-values)
    print(frequencies_table)

    scores: [] = []
    for value in frequencies_table.keys():
        # Compute Chi-Square for this value
        chi_square: float = numpy.sum(((frequencies_table[value][:] - (cycles_size * (_probabilities_xk[abs(value) - 1][:]))) ** 2) / (cycles_size * _probabilities_xk[abs(value) - 1][:]))
        # Compute the P-value for this value
        score: float = scipy.special.gammaincc(5.0 / 2.0, chi_square / 2.0)
        scores.append(score)
        
    # Return result
    elapse_time = (time.time()-start_time)*1000
    if all(score >= significance_value for score in scores):
        return ('Random Excursion', True, numpy.array(scores), elapse_time)
    return ('Random Excursion', False, numpy.array(scores), elapse_time)

def random_excursion_variant(bits: numpy.ndarray):
    """
    Random excursion variant test as described in NIST paper: https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication800-22r1a.pdf
    The focus of this test is the total number of times that a particular state is visited (i.e., occurs) in a cumulative sum random walk.
    The purpose of this test is to detect deviations from the expected number of visits to various states in the random walk.
    This test is actually a series of eighteen tests (and conclusions), one test and conclusion for each of the
    states: -9, -8, ..., -1 and +1, +2, ..., +9.

    The significance value of the test is 0.01.
    
    It is recommended that each sequence to be tested consist of a minimum of 1,000,000 bits (i.e., n ≥ 106
    """
    # Defining Test's Signficiant Value
    significance_value = 0.01
    # Defining Start Time
    start_time = time.time()
    # Copy the bits to a new array
    bits_copy: numpy.ndarray = bits.copy()
    # Convert all the zeros in the array to -1
    bits_copy[bits_copy == 0] = -1
    # Generate the padded cumulative sum of the array of -1, 1
    sum_prime: numpy.ndarray = numpy.concatenate((numpy.array([0]), numpy.cumsum(bits_copy), numpy.array([0]))).astype(int)
    # Count the number of cycles in S' (sum_prime)
    cycles_size: int = numpy.count_nonzero(sum_prime[1:] == 0)
    # Generate the counts of offsets
    unique, counts = numpy.unique(sum_prime[abs(sum_prime) < 10], return_counts=True)
    # Compute the scores (P-values)
    scores: [] = []
    for key, value in zip(unique, counts):
        # Compute the P-value for this value (if not zero)
        if key != 0:
            scores.append(abs(value - cycles_size) / math.sqrt(2.0 * cycles_size * ((4.0 * abs(key)) - 2.0)))
    # Return result
    elapse_time = (time.time()-start_time)*1000

    # Return result
    if all(score >= significance_value for score in scores):        
        return ('Random Excursion Variant', True, numpy.array(scores), elapse_time)
    return ('Random Excursion Variant', False, numpy.array(scores), elapse_time)
