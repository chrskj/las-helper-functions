import numpy as np
import matplotlib.pyplot as plt


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


def replace_values(las_object, replace_values_dict):
    for i, _ in enumerate(las_object.data):
        for curve in las_object.curves:
            for key, value in replace_values_dict.items():
                # print(f"{las_object.curves}, {curve.mnemonic}")
                # print(f"{i} {las_object[curve.mnemonic][i]}, {curve.mnemonic}")
                if key in str(las_object[curve.mnemonic][i]):
                    las_object[curve.mnemonic][i] = value

def remove_negative_values(las_object, curve_name, logging=false):
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


# Plotting ###########################################################################################################


def plot_histogram(data, labels, title, xlabel, ylabel):
    plt.hist(data, bins=20, alpha=0.5, label=labels)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()
