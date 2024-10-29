import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mode
import csv
from collections import defaultdict


def get_depth_min_max_index(las_object, depth_min, depth_max, depth_curve_name="DEPT"):
    depth_min_index = np.nan
    depth_max_index = np.nan
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
        selected_curve[i] = np.nan


def get_average(las_object, curve_name, depth_min, depth_max, depth_curve_name="DEPT"):
    depth_min_index, depth_max_index = get_depth_min_max_index(las_object, depth_min, depth_max, depth_curve_name)
    selected_curve = las_object[curve_name][depth_min_index : depth_max_index + 1]
    print(selected_curve)
    average_value = selected_curve.mean()
    return average_value


def get_averages_for_intervals(las_object, curve_name, intervals, depth_curve_name="DEPT"):
    averages = []
    for i in range(0, len(intervals) - 1):
        depth_min = intervals[i]
        depth_max = intervals[i + 1]
        avg = get_average(las_object, curve_name, depth_min, depth_max, depth_curve_name)
        averages.append((depth_min, depth_max, avg))
    return averages


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
            selected_curve[i] = np.nan


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
    return np.nan  # Handle the case where vsh is None, you might return a default value or raise an exception


def thermal_conductivity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -2.4 + (2.393 * rho_b) - (1.29 * vsh)  # A24
    return np.nan


def thermal_conductivity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return 15.69 - (3.455 * rho_b) - (0.01725 * dt)  # A7
    return np.nan


# Thermal diffusivity equations
def thermal_diffusivity_clastic(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -1.64 + (1.29 * rho_b) - (0.78 * vsh)  # B55
    return np.nan


def thermal_diffusivity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return -2.11 + (1.42 * rho_b) - (0.35 * vsh)  # B24
    return np.nan


def thermal_diffusivity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return 8.5 - (1.9 * rho_b) - (0.01065 * dt)  # B7
    return np.nan


# Specific heat capacity equations
def specific_heat_capacity_clastic(rho_b, vsh):
    if None not in (rho_b, vsh):
        return 5176.2 - (1598.4 * rho_b) - (206.8 * vsh)  # C55
    return np.nan


def specific_heat_capacity_carbonate(rho_b, vsh):
    if None not in (rho_b, vsh):
        return 5466.7 - (1664.8 * rho_b) - (436.5 * vsh)  # C24
    return np.nan


def specific_heat_capacity_evaporite(rho_b, dt):
    if None not in (rho_b, dt):
        return -1002.1 + (305 * rho_b) + (6.607 * dt)  # C7
    return np.nan


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
                thermal_conductivity = np.nan
                thermal_diffusivity = np.nan
                specific_heat_capacity = np.nan
        else:
            print(f"Missing rock type data at depth: {depth}")
            thermal_conductivity = np.nan
            thermal_diffusivity = np.nan
            specific_heat_capacity = np.nan

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

    colspans = [3, 1, 3, 3, 3, 1]
    _, axes = setup_plot(colspans)

    plot_gamma_ray(axes[0], gamma, depth)
    plot_rock_types(axes[1], well_log, rock_types, depth)
    plot_temperature(axes[2], temp, depth_discrete)
    plot_bulk_thermal_conductivity(axes[3], tc, depth_discrete)
    plot_geothermal_gradient(axes[4], gg, depth_discrete)
    plot_groups(axes[5], groups)

    finalize_plot(axes, groups, top_depth, bottom_depth)


def make_geothermal_plot_2(well_log, groups, well_log_discrete, top_depth, bottom_depth, rock_types):
    depth = well_log["DEPT"]
    depth_discrete = well_log_discrete["Depth"]
    tc = well_log_discrete["TC"]
    text_labels = [[1000, 2000, 2500], ["T = N * (T - 32) * 5/9", "T = G / (F - 32) * 5/9", "T = NR * (T - 32) * 5/9"]]

    colspans = [2, 3, 1]
    _, axes = setup_plot(colspans)

    plot_groups(axes[0], groups)
    plot_bulk_thermal_conductivity(axes[1], tc, depth_discrete, True)
    plot_rock_types(axes[2], well_log, rock_types, depth, True)
    plot_text(axes[2], text_labels)

    finalize_plot(axes, groups, top_depth, bottom_depth)


def setup_plot(colspans):
    fig, ax = plt.subplots(figsize=(sum(colspans), 10))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    columns = sum(colspans)
    rows = 1
    axes = []

    colspan_tot = 0
    for colspan in colspans:
        sharey = axes[0] if axes else None
        ax = plt.subplot2grid((rows, columns), (0, colspan_tot), rowspan=1, colspan=colspan, sharey=sharey)
        axes.append(ax)
        colspan_tot += colspan

    return fig, axes


def plot_gamma_ray(ax, gamma, depth, hide_y_labels=False):
    color_graph = "green"
    ax.plot(gamma, depth, color=color_graph, linewidth=0.5)
    ax.set_xlabel("Gamma")
    ax.set_xlim(0, 150)
    ax.set_ylabel("Depth (m)")
    ax.tick_params(axis="x", colors=color_graph)
    ax.spines["top"].set_edgecolor(color_graph)
    ax.set_xticks([0, 50, 100, 150])
    ax.set_xticklabels([])

    # Shading for Gamma-ray
    left_col_value = 0
    right_col_value = 150
    span = abs(left_col_value - right_col_value)
    cmap = plt.get_cmap("hot_r")
    color_index = np.arange(left_col_value, right_col_value, span / 100)
    for index in sorted(color_index):
        index_value = (index - left_col_value) / span
        color = cmap(index_value)
        ax.fill_betweenx(depth, gamma, right_col_value, where=gamma >= index, color=color)
    if hide_y_labels:
        ax.tick_params(labelleft=False)


def plot_rock_types(ax, well_log, rock_types, depth, hide_y_labels=False):
    # ax.xaxis.set_visible(False)
    # ax.yaxis.set_visible(False)
    color_graph = "black"
    ax.set_xlabel("Rock Types")
    ax.set_xlim(0, 1)
    ax.tick_params(axis="x", colors=color_graph)
    ax.spines["top"].set_edgecolor(color_graph)

    for key, value in rock_types.items():
        color = value["color"]
        mask = well_log["4ROCKTYPES"] == key
        ax.fill_betweenx(depth, 0, 1, where=mask, facecolor=color)

    if hide_y_labels:
        ax.tick_params(labelleft=False)
    ax.set_xticklabels([])


def plot_temperature(ax, temp, depth_discrete, hide_y_labels=False):
    color_graph = "red"
    ax.plot(temp, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax.set_xlabel("Temperature [°C]")
    ax.set_xlim(0, 200)
    ax.tick_params(axis="x", colors=color_graph)
    ax.spines["top"].set_edgecolor(color_graph)
    if hide_y_labels:
        ax.tick_params(labelleft=False)
    ax.set_xticklabels([])


def plot_bulk_thermal_conductivity(ax, tc, depth_discrete, hide_y_labels=False):
    color_graph = "purple"
    ax.plot(tc, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax.set_xlabel("Bulk thermal conduc [W/m°C]")
    ax.set_xlim(0.5, 5.5)
    ax.tick_params(axis="x", colors=color_graph)
    ax.spines["top"].set_edgecolor(color_graph)
    if hide_y_labels:
        ax.tick_params(labelleft=False)
    ax.set_xticklabels([])


def plot_geothermal_gradient(ax, gg, depth_discrete, hide_y_labels=False):
    color_graph = "brown"
    ax.plot(gg, depth_discrete, color=color_graph, linewidth=1, marker="o", markersize=2.5)
    ax.set_xlabel("Geothermal gradient [°C/km]")
    ax.set_xlim(10, 70)
    ax.tick_params(axis="x", colors=color_graph)
    ax.spines["top"].set_edgecolor(color_graph)
    if hide_y_labels:
        ax.tick_params(labelleft=False)
    ax.set_xticklabels([])


def plot_groups(ax, groups, hide_y_labels=False):
    ax.set_xlabel("Groups")
    groups_midpoints = [(value[0] + (value[1] - value[0]) / 2) for value in groups.values()]
    for label, formation_mid in zip(groups.keys(), groups_midpoints):
        ax.text(0.1, formation_mid, label, rotation=0, verticalalignment="center", fontweight="bold", fontsize="medium")
    if hide_y_labels:
        ax.tick_params(labelleft=False)
    ax.set_xticklabels([])


def plot_text(ax, text_dict):
    ax2 = ax.twinx()
    ax2.set_yticks(text_dict[0])
    ax2.set_yticklabels(text_dict[1])


def finalize_plot(axes, groups, top_depth, bottom_depth):
    for ax in axes[:-1]:  # All except the last
        ax.set_ylim(bottom_depth, top_depth)
        ax.grid(which="major", color="lightgrey", linestyle="-")
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        ax.spines["top"].set_position(("axes", 1.02))

    # Color each interval
    # for ax in axes:
    for ax in axes[:-1]:  # All except the last
        for group_i in groups.values():
            ax.axhspan(group_i[0], group_i[1], color=group_i[2], alpha=0.1)

    plt.tight_layout()
    axes[0].figure.subplots_adjust(wspace=0)
    plt.show()
