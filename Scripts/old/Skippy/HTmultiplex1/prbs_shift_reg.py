# Computes maximum length PRBS for a given polynomial
# and initial state.  The polynomial is specified by
# a list of taps.  The initial state is specified by
# a list of bits.  The function returns the PRBS

def prbs_shift_reg(taps, init_state):
    # Compute the maximum length of the shift register
    max_length = 2**len(init_state) - 1
    #print(max_length)

    # Initialize the shift register
    sr = init_state
    #print(sr)

    # Initialize the PRBS sequence
    prbs = []

    # Generate the PRBS sequence
    for i in range(max_length):
        # Compute the next state of the shift register
        next_state = 0
        #print(taps)
        for tap in taps:
            print(tap)
            print(sr[tap])
            next_state ^= sr[tap]
            #print(next_state)
        # Update the shift register
        sr = [next_state] + sr[:-1]

        # Output the next bit of the PRBS sequence
        prbs.append(sr[-1])

    # Return the PRBS sequence
    return prbs

# Test the function with the example from the prompt
taps = [2, 1]
init_state = [0, 0, 1]
prbs = prbs_shift_reg(taps, init_state)
print(prbs)
print(len(prbs))
'''
def prbs_shift_reg2(poly_coefs, init_state):
    # Compute the maximum length of the shift register
    max_length = 2**len(init_state) - 1
    #print(max_length)

    # Initialize the shift register
    sr = init_state
    #print(sr)

    # Initialize the PRBS sequence
    prbs = []

    for i in range(max_length):
        



    # Generate the PRBS sequence
    for i in range(max_length):
        # Compute the next state of the shift register
        next_state = 0
        #print(taps)
        for j, coef in enumerate(poly_coefs):
            #print(tap)
            #print(sr[tap])
            next_state ^= sr[j] * coef
            #print(next_state)
        # Update the shift register
        sr = [next_state] + sr[:-1]

        # Output the next bit of the PRBS sequence
        prbs.append(sr[-1])

    # Return the PRBS sequence
    return prbs
'''