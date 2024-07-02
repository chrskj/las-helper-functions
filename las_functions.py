import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mode
import csv
from collections import defaultdict


def get_depth_min_max_index(las_object, depth_min, depth_max, depth_curve_name="DEPT"):
    depth_min_index = np.NAN
    depth_max_index = np.NAN
    depth_curve = las_object[depth_curve_name]
    for i, depth in enumerate(depth_curve):
        if depth >= depth_min and np.isnan(depth_min_index):
            depth_min_index = i
        if depth <= depth_max:
            depth_max_index = i
        else:
            break
    return depth_min_index, depth_max_index


def remove_curve_values(las_object, curve_name, depth_min, depth_max, depth_curve_name="DEPT"):
    depth_min_index, depth_max_index = get_depth_min_max_index(las_object, depth_min, depth_max, depth_curve_name)
    selected_curve = las_object[curve_name]
    for i in range(depth_min_index, depth_max_index + 1):
        selected_curve[i] = np.NAN


def get_average(las_object, curve_name, depth_min, depth_max, depth_curve_name="DEPT"):
    depth_min_index, depth_max_index = get_depth_min_max_index(las_object, depth_min, depth_max, depth_curve_name)
    selected_curve = las_object[curve_name][depth_min_index : depth_max_index + 1]
    print(selected_curve)
    average_value = selected_curve.mean()
    return average_value


def delete_curve(las_object, curve_name):
    if curve_name in las_object.curves:
        del las_object.curves[curve_name]
        print(f"Property {curve_name} deleted.")
    else:
        print(f"Property {curve_name} doesn't exist.")


def save_las_file(las_object, filename):
    las_file_path = filename
    las_object.write(las_file_path)


def detect_step_error(las_object, depth_curve_name="DEPT"):
    step = las_object.header["Well"]["STEP"].value
    print(f"step: {step}")
    depth_curve = las_object[depth_curve_name]
    for i in range(len(depth_curve) - 1):
        next_step_val = round(depth_curve[i] + step, 4)
        print(f"{depth_curve[i + 1]} == {next_step_val}")
        if depth_curve[i + 1] != next_step_val:
            print(f"{i} - {depth_curve[i + 1]} is not correct step value")
            break


def compare_las(las_object1, las_object2):
    for curve in las_object1.curves:
        if curve in las_object2.curves:
            for i, _ in enumerate(las_object1[curve.mnemonic]):
                if las_object1[curve.mnemonic][i] != las_object2[curve.mnemonic][i] and not (
                    np.isnan(las_object1[curve.mnemonic][i]) and np.isnan(las_object2[curve.mnemonic][i])
                ):
                    print(
                        f"Curve {curve.mnemonic} is different at index {i} with values {las_object1[curve.mnemonic][i]} and {las_object2[curve.mnemonic][i]}"
                    )
        else:
            print(f"Curve {curve.mnemonic} doesn't exist in the second LAS file")

    for curve in las_object2.curves:
        if curve not in las_object1.curves:
            print(f"Curve {curve.mnemonic} doesn't exist in the first LAS file")


def replace_values(las_object, replace_values_dict):
    for i, _ in enumerate(las_object.data):
        for curve in las_object.curves:
            for key, value in replace_values_dict.items():
                # print(f"{las_object.curves}, {curve.mnemonic}")
                # print(f"{i} {las_object[curve.mnemonic][i]}, {curve.mnemonic}")
                if key in str(las_object[curve.mnemonic][i]):
                    las_object[curve.mnemonic][i] = value


def remove_negative_values(las_object, curve_name, logging=False):
    selected_curve = las_object[curve_name]
    for i, _ in enumerate(selected_curve):
        if selected_curve[i] < 0:
            if logging:
                print(f"Negative value found at index {i} with value {selected_curve[i]}")
            selected_curve[i] = np.NAN


def get_array_with_mask(las_object, array_name, mask_array_name):
    temp_array = np.array(las_object[array_name].data)
    mask_array = np.array(las_object[mask_array_name].data)
    valid_data_mask = ~(mask_array == str(las_object.well.NULL.value))
    return temp_array[valid_data_mask]


def change_null_value(las_object, new_null_value):
    current_null_value = str(las_object.header["Well"]["NULL"].value)
    las_object.header["Well"]["NULL"] = new_null_value
    replace_values_dict = {current_null_value: new_null_value}
    replace_values(las_object, replace_values_dict)


def thermal_conductivity_clastic(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -1.28 + (1.974 * rho_b) - (2.02 * vsh)  # A55
    return np.NAN  # Handle the case where vsh is None, you might return a default value or raise an exception


def thermal_conductivity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -2.4 + (2.393 * rho_b) - (1.29 * vsh)  # A24
    return np.NAN


def thermal_conductivity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return 15.69 - (3.455 * rho_b) - (0.01725 * dt)  # A7
    return np.NAN


# Thermal diffusivity equations
def thermal_diffusivity_clastic(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -1.64 + (1.29 * rho_b) - (0.78 * vsh)  # B55
    return np.NAN


def thermal_diffusivity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -2.11 + (1.42 * rho_b) - (0.35 * vsh)  # B24
    return np.NAN


def thermal_diffusivity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return 8.5 - (1.9 * rho_b) - (0.01065 * dt)  # B7
    return np.NAN


# Specific heat capacity equations
def specific_heat_capacity_clastic(rho_b, vsh):
    if None not in (rho_b, vsh):
        return 5176.2 - (1598.4 * rho_b) - (206.8 * vsh)  # C55
    return np.NAN


def specific_heat_capacity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return 5466.7 - (1664.8 * rho_b) - (436.5 * vsh)  # C24
    return np.NAN


def specific_heat_capacity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return -1002.1 + (305 * rho_b) + (6.607 * dt)  # C7
    return np.NAN


def calculate_thermal_properties(rho_b, dt, vsh, rock_type, rock_types):
    if rock_types[rock_type]["rock_type"] in ("Sandstones", "Claystones"):
        return (
            thermal_conductivity_clastic(rho_b, vsh),
            thermal_diffusivity_clastic(rho_b, vsh),
            specific_heat_capacity_clastic(rho_b, vsh),
        )
    if rock_types[rock_type]["rock_type"] == "Carbonates":
        return (
            thermal_conductivity_carbonate(rho_b, vsh),
            thermal_diffusivity_carbonate(rho_b, vsh),
            specific_heat_capacity_carbonate(rho_b, vsh),
        )
    if rock_types[rock_type]["rock_type"] == "Evaporites":
        return (
            thermal_conductivity_evaporite(rho_b, dt),
            thermal_diffusivity_evaporite(rho_b, dt),
            specific_heat_capacity_evaporite(rho_b, dt),
        )
    raise ValueError(f"Unknown rock_type: {rock_type}")


def add_thermal_properties(
    las_object,
    rock_types,
    curve_names,
    depth_curve_name="DEPT",
):
    # Extract well log data
    thermal_conductivity_data = []
    thermal_diffusivity_data = []
    specific_heat_capacity_data = []
    # Calculate TC, TD, and SHC for each depth point
    for i, depth in enumerate(las_object[depth_curve_name]):
        rock_type = las_object[curve_names["rock_type"]][i]
        if not np.isnan(rock_type):
            rho_b = las_object[curve_names["rho_b"]][i]
            phi_n = las_object[curve_names["phi_n"]][i]
            dt = las_object[curve_names["dt"]][i]
            vsh = las_object[curve_names["vsh"]][i]
            # print(f"{depth}: {rho_b} {phi_n} {dt} {vsh}{lithology}")

            if all((rho_b, dt, vsh)):
                thermal_conductivity, thermal_diffusivity, specific_heat_capacity = calculate_thermal_properties(
                    rho_b, dt, vsh, rock_type, rock_types
                )
                # print(thermal_conductivity)
            else:
                print(f"Missing well log data at depth: {depth} {rho_b}, {phi_n}, {dt}, {vsh}")
                thermal_conductivity = np.NAN
                thermal_diffusivity = np.NAN
                specific_heat_capacity = np.NAN
        else:
            print(f"Missing rock type data at depth: {depth}")
            thermal_conductivity = np.NAN
            thermal_diffusivity = np.NAN
            specific_heat_capacity = np.NAN

        # print(f"{thermal_conductivity} {thermal_diffusivity} {specific_heat_capacity}")
        thermal_conductivity_data.append(thermal_conductivity)
        thermal_diffusivity_data.append(thermal_diffusivity)
        specific_heat_capacity_data.append(specific_heat_capacity)

    # Append the calculated values as new curves to the existing LAS file
    las_object.append_curve("TC", unit="W/mK", data=thermal_conductivity_data)
    las_object.append_curve("TD", unit="m2/s", data=thermal_diffusivity_data)
    las_object.append_curve("SHC", unit="J/kgK", data=thermal_diffusivity_data)
    # las.append_curve(f"SHC", unit="J/kgK", data=specific_heat_capacity_data)


def load_groups(groups_file_path):
    group_dict = {}
    with open(groups_file_path, encoding="utf-8") as file:
        next(file)  # skip header row
        for row in csv.DictReader(file, fieldnames=["Group", "Top", "Bottom", "Zone Colors"], delimiter=";"):
            group_dict[row["Group"]] = [float(row["Top"]), float(row["Bottom"]), row["Zone Colors"]]
    return group_dict


# Plotting ###########################################################################################################


def plot_histogram(data, labels, title, xlabel, ylabel):
    plt.hist(data, bins=20, alpha=0.5, label=labels)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()


def make_geothermal_plot(well_log, groups, well_log_discrete, top_depth, bottom_depth, rock_types):
    depth = well_log["DEPT"]
    gamma = well_log["HGR(NEW)"]
    depth_discrete = well_log_discrete["Depth"]
    tc = well_log_discrete["TC"]
    gg = well_log_discrete["GG"]
    temp = well_log_discrete["Temp"]

    fig, ax = plt.subplots(figsize=(15, 10))
    # Remove unused tick labels for whole graph
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Set up the plot axes
    rows, columns = 1, 14
    rows_i = 0
    colspan_i = 3
    colspan_tot = 0
    ax1 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i)
    colspan_tot += colspan_i
    colspan_i = 1
    ax2 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i, sharey=ax1)
    colspan_tot += colspan_i
    colspan_i = 3
    ax3 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i, sharey=ax1)
    colspan_tot += colspan_i
    colspan_i = 3
    ax4 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i, sharey=ax1)
    colspan_tot += colspan_i
    colspan_i = 3
    ax5 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i, sharey=ax1)
    colspan_tot += colspan_i
    colspan_i = 1
    ax6 = plt.subplot2grid((rows, columns), (rows_i, colspan_tot), rowspan=1, colspan=colspan_i, sharey=ax1)

    # As our curve scales will be detached from the top of the track,
    # this code adds the top border back in without dealing with splines
    ax10 = ax1.twiny()
    ax10.xaxis.set_visible(False)
    ax11 = ax2.twiny()
    ax11.xaxis.set_visible(False)
    ax12 = ax3.twiny()
    ax12.xaxis.set_visible(False)
    ax13 = ax4.twiny()
    ax13.xaxis.set_visible(False)
    ax14 = ax5.twiny()
    ax14.xaxis.set_visible(False)

    # Gamma-ray track ========================================================

    color_graph = "green"
    ax1.plot(gamma, depth, color=color_graph, linewidth=0.5)
    ax1.set_xlabel("Gamma")
    ax1.xaxis.label.set_color(color_graph)
    ax1.set_xlim(0, 150)
    ax1.set_ylabel("Depth (m)")
    ax1.tick_params(axis="x", colors=color_graph)
    ax1.spines["top"].set_edgecolor(color_graph)
    ax1.title.set_color(color_graph)
    ax1.set_xticks([0, 50, 100, 150])
    ax1.text(0.05, 1.04, 0, color=color_graph, horizontalalignment="left", transform=ax1.transAxes)
    ax1.text(0.95, 1.04, 150, color=color_graph, horizontalalignment="right", transform=ax1.transAxes)
    ax1.set_xticklabels([])

    ## Setting Up Shading for GR
    left_col_value = 0
    right_col_value = 150
    span = abs(left_col_value - right_col_value)
    cmap = plt.get_cmap("hot_r")
    color_index = np.arange(left_col_value, right_col_value, span / 100)
    # loop through each value in the color_index
    for index in sorted(color_index):
        index_value = (index - left_col_value) / span
        color = cmap(index_value)  # obtain color for color index value
        ax1.fill_betweenx(depth, gamma, right_col_value, where=gamma >= index, color=color)

    # Rock type track ========================================================

    ax2.set_xlabel("Rock Types")
    ax2.set_xlim(0, 1)
    ax2.xaxis.label.set_color("black")
    ax2.tick_params(axis="x", colors="black")
    ax2.spines["top"].set_edgecolor("black")

    for key, value in rock_types.items():
        color = value["color"]
        mask = well_log["4ROCKTYPES"] == key
        ax2.fill_betweenx(depth, 0, 1, where=mask, facecolor=color)
    ax2.set_xticklabels([])

    # Temperature track ========================================================

    color_graph = "red"
    ax3.plot(temp, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax3.set_xlabel("Temperature [°C]")
    ax3.set_xlim(0, 200)
    ax3.xaxis.label.set_color(color_graph)
    ax3.tick_params(axis="x", colors=color_graph)
    ax3.spines["top"].set_edgecolor(color_graph)
    ax3.text(0.05, 1.04, 0, color=color_graph, horizontalalignment="left", transform=ax3.transAxes)
    ax3.text(0.95, 1.04, 200, color=color_graph, horizontalalignment="right", transform=ax3.transAxes)
    ax3.set_xticklabels([])

    # Bulk thermal conductivity ========================================================

    color_graph = "purple"
    ax4.plot(tc, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax4.set_xlabel("Bulk thermal conduc [W/m°C]")
    ax4.set_xlim(0.5, 5.5)
    ax4.xaxis.label.set_color(color_graph)
    ax4.tick_params(axis="x", colors=color_graph)
    ax4.spines["top"].set_edgecolor(color_graph)
    ax4.text(0.05, 1.04, 0, color=color_graph, horizontalalignment="left", transform=ax4.transAxes)
    ax4.text(0.95, 1.04, 5.5, color=color_graph, horizontalalignment="right", transform=ax4.transAxes)
    ax4.set_xticklabels([])

    # Geothermal gradient ========================================================

    color_graph = "brown"
    ax5.plot(gg, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax5.set_xlabel("Geothermal gradient [°C/km]")
    ax5.set_xlim(10, 70)
    ax5.xaxis.label.set_color(color_graph)
    ax5.tick_params(axis="x", colors=color_graph)
    ax5.spines["top"].set_edgecolor(color_graph)
    ax5.text(0.05, 1.04, 0, color=color_graph, horizontalalignment="left", transform=ax5.transAxes)
    ax5.text(0.95, 1.04, 70, color=color_graph, horizontalalignment="right", transform=ax5.transAxes)
    ax5.set_xticklabels([])

    # Formations ========================================================

    ax6.set_xticklabels([])
    ax6.text(0.5, 1.1, "Formations", fontweight="bold", horizontalalignment="center", transform=ax6.transAxes)

    # Adding in neutron density shading

    # Common functions for setting up the plot can be extracted into
    # a for loop. This saves repeating code.
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_ylim(bottom_depth, top_depth)
        ax.grid(which="major", color="lightgrey", linestyle="-")
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        ax.spines["top"].set_position(("axes", 1.02))

    for ax in [ax1, ax3, ax4, ax5, ax6]:
        # loop through the formations dictionary and zone colors
        for group_i in groups.values():
            # use the depths and colors to shade across the subplots
            ax.axhspan(group_i[0], group_i[1], color=group_i[2], alpha=0.1)

    for ax in [ax2, ax3, ax4, ax5, ax6]:
        plt.setp(ax.get_yticklabels(), visible=False)

    groups_midpoints = []
    for _, value in groups.items():
        groups_midpoints.append(value[0] + (value[1] - value[0]) / 2)

    for label, formation_mid in zip(groups.keys(), groups_midpoints):
        ax6.text(
            0.2, formation_mid, label, rotation=45, verticalalignment="center", fontweight="bold", fontsize="medium"
        )

    plt.tight_layout()
    fig.subplots_adjust(wspace=0)
    plt.show()
