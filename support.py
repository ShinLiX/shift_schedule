# support.py
"""
Support functions for the employee scheduling system.
"""

def compute_total_weekly_availability(employees_availability):
    """
    Return a dictionary mapping each employee to the total
    number of available weekly slots (0..83).
    """
    return {e: len(avail) for e, avail in employees_availability.items()}


def is_contiguous_ones(bits):
    """
    Check if a list of bits consists of exactly one contiguous
    block of 1s or is all zero.
    """
    idxs = [i for i, b in enumerate(bits) if b == 1]
    if not idxs:
        return True
    first = idxs[0]
    last = idxs[-1]
    return (last - first + 1) == len(idxs)


def generate_day_patterns(availability_mask, min_block=6, max_block=8):
    """
    Generate valid single-day patterns (12 hours => indices 0..11).
    
    Logic:
      - sum=0 => off
      - if total availability < 6 => sum=0 or exactly fill the smaller block (contiguous)
      - else => one contiguous block in [6..8], or off

    Returns a list of patterns, each pattern is length=13 => 
    [bit0..bit11, workedBit], where workedBit=1 if any bit=1, else 0.
    """
    n = sum(availability_mask)
    valid_patterns = []

    # Iterate over all possible 12-bit patterns (2^12)
    for mask in range(1 << 12):
        pat = [(mask >> i) & 1 for i in range(12)]

        # Must be subset of availability
        if any(pat[i] == 1 and availability_mask[i] == 0 for i in range(12)):
            continue

        total_ones = sum(pat)
        worked_bit = 1 if total_ones > 0 else 0

        if n < min_block:
            # If total availability < 6 => sum=0 or fill entire smaller block if contiguous
            if total_ones == 0:
                valid_patterns.append(pat + [0])
            elif total_ones == n:
                # Must match availability exactly & be contiguous
                if pat == availability_mask and is_contiguous_ones(pat):
                    valid_patterns.append(pat + [1])
        else:
            # availability >= 6 => sum=0 => off, or single contiguous block in [6..8]
            if total_ones == 0:
                valid_patterns.append(pat + [0])
            else:
                if (min_block <= total_ones <= max_block) and is_contiguous_ones(pat):
                    valid_patterns.append(pat + [1])

    return valid_patterns