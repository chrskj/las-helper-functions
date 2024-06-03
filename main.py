import numpy as np
import lasio
from las_functions import *

lithologies = {
    0: {"lith": "Fine Sand", "lith_num": 0, "hatch": "...", "color": "#ffffbf"},
    1: {"lith": "Sandstone", "lith_num": 1, "hatch": "..", "color": "#ffff00"},
    2: {"lith": "Sandstone/Shale", "lith_num": 2, "hatch": "-.", "color": "#ffe119"},
    3: {"lith": "Shale", "lith_num": 3, "hatch": "--", "color": "#bebebe"},
    4: {"lith": "Marl", "lith_num": 4, "hatch": "", "color": "#7cfc00"},
    5: {"lith": "Dolomite", "lith_num": 5, "hatch": "-/", "color": "#8080ff"},
    6: {"lith": "Limestone", "lith_num": 6, "hatch": "+", "color": "#80ffff"},
    7: {"lith": "Chalk", "lith_num": 7, "hatch": "..", "color": "#80ffff"},
    8: {"lith": "Halite", "lith_num": 8, "hatch": "x", "color": "#7ddfbe"},
    9: {"lith": "Anhydrite", "lith_num": 9, "hatch": "", "color": "#ff80ff"},
    10: {"lith": "Tuff", "lith_num": 10, "hatch": "||", "color": "#ff8c00"},
    11: {"lith": "Coal", "lith_num": 11, "hatch": "", "color": "black"},
    12: {"lith": "Basement", "lith_num": 12, "hatch": "-|", "color": "#ef138a"},
}

rock_types = {
    0: {"rock_type": "Sandstones", "rock_type_num": 0, "hatch": "...", "color": "#ffffbf"},
    1: {"rock_type": "Claystones", "rock_type_num": 1, "hatch": "..", "color": "#ffff00"},
    2: {"rock_type": "Carbonates", "rock_type_num": 2, "hatch": "-.", "color": "#ffe119"},
    3: {"rock_type": "Evaporites", "rock_type_num": 3, "hatch": "--", "color": "#bebebe"},
}

LOG_SOURCE_FOLDER = "log_files"
LOG_OUTPUT_FOLDER = "log_files_output"


def get_rock_type_by_name(rock_type_name):
    for _, value in rock_types.items():
        if value["rock_type"] == rock_type_name:
            return value
    return None


# Example of
# - how to detect step errors in a given curve and print the depth values
# - how to remove curve values in a given range 
# - how to delete a curve
# - how to change null values
def main1():
    file_name = "1_3-2_Lith_Vsh_TP_RHP"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    las = lasio.read(las_file_path)
    detect_step_error(las)
    remove_curve_values(las, "HRHOB(NEW)", 1804.3124, 2093.1104)
    delete_curve(las, "LITHOLOGY")
    change_null_value(las, "-999.25000")
    save_las_file(las, f"{LOG_OUTPUT_FOLDER}/{file_name}_modified.las")

# Example of
# - how to replace values in a curve
def main2():
    file_name = r"1_3-2_Lith_Vsh_TP_RHP"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    replace_values_dict = {
        "sandstones": 0,
        "claystones": 1,
        "carbonates": 2,
        "evaporites": 3,
    }
    las = lasio.read(las_file_path)
    replace_values(las, replace_values_dict)
    save_las_file(las, f"{LOG_OUTPUT_FOLDER}/{file_name}_modified.las")

# Example of
# - how to get the average of a curve in a given range
def main3():
    file_name = r"1_2-2"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    las = lasio.read(las_file_path)
    average = get_average(las, "DEN(NEW)", 1804.17, 1807.0656)
    print(f"average is {average}")

# Example of
# - how to add new curves with thermal properties to a LAS file
def main4():
    file_name = r"2_2-1_New"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    las = lasio.read(las_file_path)
    curve_names = {
        "rock_type": "ROCKTYPES",
        "rho_b": "HRHOB(NEW)",
        "phi_n": "HNPHI(NEW)",
        "dt": "HDT(NEW)",
        "vsh": "VSH",
    }
    add_thermal_properties(las, rock_types, curve_names)
    save_las_file(las, f"{LOG_OUTPUT_FOLDER}/{file_name}_modified.las")

# Example of
# - how to plot histograms of thermal conductivity values for different lithologies
def main5():
    file_name = "1_3-2_Lith_Vsh_Thermal_Prop"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    las = lasio.read(las_file_path)

    # Accessing data from LAS file

    tc = get_array_with_mask(las, "TC", "TC")
    rock_types_vals = get_array_with_mask(las, "ROCKTYPES", "TC")

    # Grouping thermal conductivity values by lithology
    claystones = tc[rock_types_vals == get_rock_type_by_name("Claystones")["rock_type_num"]]
    sandstones = tc[rock_types_vals == get_rock_type_by_name("Sandstones")["rock_type_num"]]
    carbonates = tc[rock_types_vals == get_rock_type_by_name("Carbonates")["rock_type_num"]]
    evaporites = tc[rock_types_vals == get_rock_type_by_name("Evaporites")["rock_type_num"]]

    # Plotting histograms
    plot_histogram(
        [claystones, sandstones, carbonates, evaporites],
        labels=["Claystones", "Sandstones", "Carbonates", "Evaporites"],
        title="Thermal Conductivity Histogram",
        xlabel="Thermal Conductivity (W/m°K)",
        ylabel="Frequency",
    )

def main6():
    file_name = "1_3-1"
    las_file_path = f"{LOG_SOURCE_FOLDER}/{file_name}.las"
    las = lasio.read(las_file_path)

    remove_negative_values(las, "TD", logging=true)

    save_las_file(las, f"{LOG_OUTPUT_FOLDER}/{file_name}_modified.las")

if __name__ == "__main__":
    main6()
