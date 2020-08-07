## Pseudo Code for Wind Erosion Model
## Written by Charlie Weil & Isita Talukdar, August 2020
## Inspired by Huang Binbin's approach.

def execute(args):
	## See typical InVEST Model structure.
	return None

# Step 1 : Weather Factor 
# # # # # # # # # # # # #

def calculate_air_density(temperature, air_pressure):
	"""Calculate air density for a given month.

    Parameters:
        Temperature for a given month, in [UNIT °C??]
        Air

    Returns:
        Air density (rho), in UNIT !!
    """
	air_density_rho = 1.293 * (273 / (273 + temperature_month)) * air_pressure / 101.3

    return air_density_rho

def calculate_monthly_air_density(dem):
	"""Calculate air density for each month

    Parameters:
        Temperature for all months, in [UNIT °C??]
        DEM

    Returns:
        Air density (rho), in UNIT !!
    """
    