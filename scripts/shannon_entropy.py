#########################################################
# Plots the Shannon entropy of the available statepoint #
#########################################################

import openmc
import matplotlib.pyplot as plt

sp = openmc.StatePoint('statepoint.100.h5')
entropy = sp.entropy

# Plot the Shannon entropy
plt.plot(entropy)
plt.xlabel('Batch')
plt.ylabel('Shannon Entropy')
plt.savefig('shannon_entropy.png', format='png', dpi=300)
