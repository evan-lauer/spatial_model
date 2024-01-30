import pandas as pd
import numpy as np
import bisect
import warnings
# We have some annoying warnings that come up, so I'm suppressing them.
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)


fungi_names = pd.read_csv('./fungal_biogeography/fungi_data/Fungal_climate_data.csv', index_col=0).loc[:,['gen.name2']]

fungi_temperature_curves = pd.read_csv('./fungal_biogeography/fungi_data/Fungi_temperature_curves.csv', index_col=0)
fungi_temperature_curves = fungi_temperature_curves.loc[fungi_temperature_curves['type']=='smoothed']

fungi_moisture_curves = pd.read_csv('./fungal_biogeography/fungi_data/Fungi_moisture_curves.csv', index_col=0)
fungi_moisture_curves = fungi_moisture_curves.loc[fungi_moisture_curves['type']=='smoothed']

fungi_decomposition_chart = pd.read_csv('./fungal_biogeography/fungi_data/Fungal_decomposition_data_converted.csv', index_col=0)

max_extension_rates = {}
for index, row in fungi_names.iterrows():
  single_fungus_moisture_curve = fungi_moisture_curves.loc[[index]]
  extension_rates = single_fungus_moisture_curve.loc[[index],['hyphal_rate']]['hyphal_rate'].tolist()
  max_extension_rates[index] = max(extension_rates)

# Return the max HER divided by the actual HER by moisture
def get_moisture_strength_factor(shortname, water_potential):
  strength_factor = get_hyphal_ext_rate_by_moisture(shortname, water_potential) / max_extension_rates[shortname]
  return max(strength_factor, .1)


def get_decomposition_rate(shortname, temperature):
  if temperature == "low":
    return fungi_decomposition_chart.loc[shortname,"low_temp"]
  if temperature == "medium":
    return fungi_decomposition_chart.loc[shortname,"medium_temp"]
  if temperature == "high":
    return fungi_decomposition_chart.loc[shortname,"high_temp"]
  if temperature == "mean":
    return fungi_decomposition_chart.loc[shortname,"mean"]
  raise ValueError

  

def get_fungus_name(shortname):
  return fungi_names.loc[[shortname]].values[0][0]

# This maps the given <val> to the nearest value in <values>
# Assumes <values> is sorted in <order>
def get_nearest_value(values: list, val, order):
  if val in values:
    return val
  if order == "ascending":
    i = bisect.bisect_left(values, val)
    if i:
      return values[i-1]
    raise ValueError
  if order == "descending":
    values.reverse()
    i = bisect.bisect_left(values, val)
    if i:
      closest_val = values[i-1]
      values.reverse() # We need to un-reverse the list before we return 
      return closest_val
  raise ValueError("Order must be ascending or descending")

# Will get the hyphal ext rate of the shortnamed fungus at the temperature CLOSEST
# to <temperature>. There may not be data for the given temperature.
def get_hyphal_ext_rate_by_temperature(shortname, temperature):
  single_fungus_temperature_curve = fungi_temperature_curves.loc[[shortname]]
  indices = single_fungus_temperature_curve.loc[[shortname],['temp_c']]['temp_c'].tolist()
  index = get_nearest_value(indices, temperature, "ascending")
  row = single_fungus_temperature_curve.loc[single_fungus_temperature_curve['temp_c']==index]
  hyphal_ext_rate = row['hyphal_rate'][0]
  return float(hyphal_ext_rate)

# Will get the hyphal ext rate of the shortnamed fungus at the water potential CLOSEST
# to <water_potential>.
# Data exists for water potential of -5 to 0.
def get_hyphal_ext_rate_by_moisture(shortname, water_potential):
  single_fungus_moisture_curve = fungi_moisture_curves.loc[[shortname]]
  indices = single_fungus_moisture_curve.loc[[shortname],['matric_pot']]['matric_pot'].tolist()
  index = get_nearest_value(indices, water_potential, "descending")
  row = single_fungus_moisture_curve.loc[single_fungus_moisture_curve['matric_pot']==index]
  hyphal_ext_rate = row['hyphal_rate'][0]
  return float(hyphal_ext_rate)



# This code scales the decomposition rate over 122 days into the proper per-day decomposition rate.
# We only needed to run it once and store the results.

# def convert_decomposition_rate(rate):
#   proportion = 1 - (rate / 100)
#   return 1 - proportion ** (1 / 122)

# fungi_decomposition_chart = pd.read_csv('./fungal_biogeography/fungi_data/Fungal_decomposition_data.csv', index_col=0)
# for index, row in fungi_decomposition_chart.iterrows():
#   fungi_decomposition_chart.loc[index, "low_temp"] = convert_decomposition_rate(row["low_temp"])
#   fungi_decomposition_chart.loc[index, "medium_temp"] = convert_decomposition_rate(row["medium_temp"])
#   fungi_decomposition_chart.loc[index, "high_temp"] = convert_decomposition_rate(row["high_temp"])
#   fungi_decomposition_chart.loc[index, "mean"] = convert_decomposition_rate(row["mean"])

# fungi_decomposition_chart.to_csv('./fungal_biogeography/fungi_data/Fungal_decomposition_data_converted.csv')