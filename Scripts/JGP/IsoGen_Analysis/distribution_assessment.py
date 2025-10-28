import numpy as np
import matplotlib.pyplot as plt

def calculate_pearson_chisquared(predicted_dist, expected_dist, cutoff=0.005):
    """
    Calculate the Pearson chi-squared statistic for two distributions.

    :param predicted_dist: The predicted distribution (numpy array).
    :param expected_dist: The expected distribution (numpy array).
    :return: The Pearson chi-squared statistic.
    """

    # Ensure both distributions are the same length
    if len(predicted_dist) != len(expected_dist):
        #Pad the shorter distribution with zeros
        if len(predicted_dist) < len(expected_dist):
            predicted_dist = np.pad(predicted_dist, (0, len(expected_dist) - len(predicted_dist)))
        else:
            expected_dist = np.pad(expected_dist, (0, len(predicted_dist) - len(expected_dist)))

    # Ensure both distributions are numpy arrays
    predicted_dist = np.array(predicted_dist)
    expected_dist = np.array(expected_dist)

    #Normalize the distributions to max 1
    predicted_dist = predicted_dist / np.max(predicted_dist)
    expected_dist = expected_dist / np.max(expected_dist)

    # Apply the cutoff to both distributions
    predicted_dist[predicted_dist < cutoff] = 0
    expected_dist[expected_dist < cutoff] = 0

    chi_squared = 0

    for i in range(len(expected_dist)):
        # Skip if the expected value is zero to avoid division by zero
        if expected_dist[i] == 0:
            continue

        chi_squared += ((predicted_dist[i] - expected_dist[i])**2) / expected_dist[i]

    if False:
        plt.plot(predicted_dist, label='Predicted')
        plt.plot(expected_dist, label='Expected')
        plt.legend()
        plt.title(f'Chi-squared: {chi_squared_manual:.2f}')
        plt.show()

    return chi_squared

def calculate_cosinesimilarity(predicted_dist, expected_dist):
    """
    Calculate the cosine similarity between two distributions.

    :param predicted_dist: The predicted distribution (numpy array).
    :param expected_dist: The expected distribution (numpy array).
    :return: The cosine similarity.
    """

    # Ensure both distributions are the same length
    # Ensure both distributions are the same length
    if len(predicted_dist) != len(expected_dist):
        # Pad the shorter distribution with zeros
        if len(predicted_dist) < len(expected_dist):
            predicted_dist = np.pad(predicted_dist, (0, len(expected_dist) - len(predicted_dist)))
        else:
            expected_dist = np.pad(expected_dist, (0, len(predicted_dist) - len(expected_dist)))

    # Ensure both distributions are numpy arrays
    predicted_dist = np.array(predicted_dist)
    expected_dist = np.array(expected_dist)

    # Normalize the distributions to max 1
    predicted_dist = predicted_dist / np.max(predicted_dist)
    expected_dist = expected_dist / np.max(expected_dist)

    # Calculate the dot product of the two distributions
    dot_product = np.dot(predicted_dist, expected_dist)

    # Calculate the magnitude of each distribution
    magnitude_predicted = np.linalg.norm(predicted_dist)
    magnitude_expected = np.linalg.norm(expected_dist)

    # Calculate the cosine similarity
    cosine_similarity = dot_product / (magnitude_predicted * magnitude_expected)

    return cosine_similarity



