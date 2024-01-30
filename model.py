import pandas
import numpy as np
import matplotlib.pyplot as plt
import math
import warnings
import random
from data_manipulation import get_hyphal_ext_rate_by_moisture, get_fungus_name, get_decomposition_rate, get_moisture_strength_factor
from moisture import get_moisture_array
# We have some annoying warnings that come up, so I'm suppressing them.
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)

WATER_POTENTIAL_ARRAY = get_moisture_array()


# Parameters for the Euler approximation of the dif eqs
# Ideally SIMULATION_DURATION is divisible by STEP_SIZE
STEP_SIZE = 1                                     # step size: days
SIMULATION_DURATION = 365 * 2                         # simulation length: days
NUM_STEPS = int(SIMULATION_DURATION / STEP_SIZE)  # number of steps
NUM_FUNGI = 4                                     # number of fungi we are modeling
FUNGI_NAMES = ["p.rufa.acer.n", "p.gilv.n", "p.har.n", "p.robin.n"]       # names of the fungi we are modeling
CLIMATE_TYPE = "medium"                           # either "low", "medium", "high"

# Fixed constants to be used as parameters within our dif eqs
SURFACE_INITIAL_COVER = 500000             # starting amount of ground cover: mm^2
SURFACE_SIZE = 1000000                     # simulation surface size: mm^2
SURFACE_NEW_COVER = 200                    # new ground cover of dead matter per day: mm^2/day
FUNGI_INITIAL_RADIUS = 1                   # starting radius of each fungus: mm
INITIAL_WATER_POTENTIAL = -1.0             # starting water potential of the system: MPa

# These are variables we're modeling (ie: these are the variables that we 
# have dif eqs for)
ground_cover = SURFACE_INITIAL_COVER                               # organic ground cover: mm^2
fungal_radii = np.full(NUM_FUNGI, FUNGI_INITIAL_RADIUS, dtype="float")    # list of hyphal radii for each fungus: mm
water_potential = INITIAL_WATER_POTENTIAL                          # water potential of the system: MPa


# This initializes our plot array. We need to keep track of every fungus' radius at every time step,
# so we use a 2-dimensional array. There's an array for every time step of length NUM_FUNGI.
fungal_radii_historical = np.ndarray(shape=(NUM_STEPS, NUM_FUNGI))     # 2d array, where each array[i] represents every fungus' radii at time step i

ground_cover_historical = np.zeros(SIMULATION_DURATION)

# This initializes our decomposition array. For every fungus, it represents the total quantity
# of wood that the fungus could possibly consume in a day
DECOMPOSITION_MODE = "mean"
fungal_decomposition = np.zeros(NUM_FUNGI)    # for each fungi, max quantity of wood that could be consumed
for i in range(NUM_FUNGI):
  fungal_decomposition[i] = get_decomposition_rate(FUNGI_NAMES[i], DECOMPOSITION_MODE)
total_fungal_decomposition = sum(fungal_decomposition)

# TODO: Create the full size combat matrix (ie, (i,j) represents the outcome of i fighting j)
# for now we use a temp matrix
combat_matrix = np.array([[0, 1, 1, 1],[-1,-1, 0, -1],[-1, 0, 1, 0],[-1, -1, 1, -1]]) # this represents the fact the p.rufa beats p.pend in direct combat trials

# This is the main Euler approx. loop

for day in range(SIMULATION_DURATION):
  ground_cover_historical[day] = ground_cover
  water_potential = WATER_POTENTIAL_ARRAY[day % 365]
  decomposition_term = 0
  for i in range(NUM_FUNGI):
    # This calculates d_i for each mushroom (the product of the ground cover, the proportion of decomposition potential, and the
    # weighted moisture strength term)
    decomposition_term += (ground_cover / SURFACE_SIZE) * get_moisture_strength_factor(FUNGI_NAMES[i], water_potential)
  d_ground_cover = SURFACE_NEW_COVER - decomposition_term
  d_fungal_radii = np.zeros(NUM_FUNGI)
  for i in range(NUM_FUNGI):
    if fungal_radii[i] <= 0:
      fungal_radii[i] = FUNGI_INITIAL_RADIUS
    d_fungal_radii[i] = get_hyphal_ext_rate_by_moisture(FUNGI_NAMES[i], water_potential) * get_moisture_strength_factor(FUNGI_NAMES[i], water_potential)
    for j in range(NUM_FUNGI):
      if i != j:
        overlap_factor = 0


        if combat_matrix[i][j] == 0:
          if (random.randint(1,10) < 2):
            overlap_factor = .015 * math.pi * (fungal_radii[j] ** 2) * (fungal_radii[i] ** 2) / (SURFACE_SIZE)
          else:
            overlap_factor = .00525 *math.pi * (fungal_radii[j] ** 2) * (fungal_radii[i] ** 2) / (SURFACE_SIZE)
        elif combat_matrix[i][j] == -1:
          if (random.randint(1,10) < 2):
            overlap_factor = .00525 *math.pi * (fungal_radii[j] ** 2) * (fungal_radii[i] ** 2) / (SURFACE_SIZE)
          else:
            overlap_factor = .015 * math.pi * (fungal_radii[j] ** 2) * (fungal_radii[i] ** 2) / (SURFACE_SIZE)


        d_fungal_radii[i] -= overlap_factor * get_hyphal_ext_rate_by_moisture(FUNGI_NAMES[i], water_potential) * get_moisture_strength_factor(FUNGI_NAMES[i], water_potential) 
    
  fungal_radii_historical[day][:] = fungal_radii                           # set the first historical entry for time step 0

  ground_cover += d_ground_cover
  np.add(fungal_radii, d_fungal_radii, out=fungal_radii,casting="unsafe")
  for i in range(NUM_FUNGI):
    if fungal_radii[i] < 0:
      fungal_radii[i] = FUNGI_INITIAL_RADIUS
    elif fungal_radii[i] > 564:
      fungal_radii[i] = 560

  

plottable = np.transpose(fungal_radii_historical)
plt.subplot(2, 1, 1)  # 2 rows, 1 column, first subplot


plt.plot(plottable[0], label='Fungus 1')
plt.plot(plottable[0] + plottable[1], label='Fungus 2')
plt.plot(plottable[0] + plottable[1] + plottable[2], label='Fungus 3')
plt.plot(plottable[0] + plottable[1] + plottable[2] + plottable[3], label='Fungus 4')

plt.fill_between(range(len(plottable[0])), 0, plottable[0], alpha=0.3, color='blue')
plt.fill_between(range(len(plottable[0])), plottable[0], plottable[0] + plottable[1], alpha=0.3, color='orange')
plt.fill_between(range(len(plottable[0])), plottable[0] + plottable[1], plottable[0] + plottable[1] + plottable[2], alpha=0.3, color='green')
plt.fill_between(range(len(plottable[0])), plottable[0] + plottable[1] + plottable[2], plottable[0] + plottable[1] + plottable[2] + plottable[3], alpha=0.3, color='red')

plt.xlabel('Days Elapsed')
plt.ylabel('Hyphal Radius (mm)')
plt.title("Hyphal Extension Rate Over Time")


plt.legend()

plt.subplot(2, 1, 2)  # 2 rows, 1 column, second subplot
plt.plot(ground_cover_historical, label='Ground Cover')
plt.xlabel('Days Elapsed')
plt.ylabel('Ground Cover (mm^2)')

plt.show()