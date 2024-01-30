import numpy as np

def get_moisture_array():
  # Winter: December to February
  winter_values = np.random.uniform(low=-1.0, high=0.0, size=90)

  # Spring: March to May
  spring_values = np.random.uniform(low=-2.9, high=-.7, size=92)

  # Summer: June to August
  summer_values = np.random.uniform(low=-1.2, high=0.0, size=92)

  # Fall: September to November
  fall_values = np.random.uniform(low=-3.1, high=-1.7, size=91)

  return np.concatenate([winter_values, spring_values, summer_values, fall_values])
