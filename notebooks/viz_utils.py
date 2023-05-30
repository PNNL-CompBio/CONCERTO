import cobra
import os
import json
import math
import pandas as pd
import numpy as np
import d3flux as d3f


# Function to build combined model by looping through external metabolites
# Note:
#   - corresponding reactions named with a "EX_"+metab_id+org_id format
#   -


def make_combined_model_external_mets(org_list, model_name):
    metabolite_exchange_model = cobra.Model(model_name)
    # go through each organism
    for this_org in org_list:
        with this_org as org:
            # create new metabolite species (needed for d3flux)
            org_m = cobra.Metabolite(org.id)
            org_m.compartment = org.id

            # go through each non-external metabolite
            org_m_reactions = []
            org_em = [m for m in org.metabolites if m.compartment == 'e']
            for m in org_em:
                # Add new reaction for this metabolite
                #                 new_r = cobra.Reaction("Reaction_" + m.id + "_" + org.id)
                new_r = cobra.Reaction("EX_" + m.id + "_" + org.id)
                new_r.add_metabolites({org_m: -1, m: 1})
                org_m_reactions.append(new_r)

            # Add new reactions to new model
            metabolite_exchange_model.add_reactions(org_m_reactions)

    return metabolite_exchange_model


def make_combined_model_external_mets_shared_only(org_list, model_name):
    metabolite_exchange_model = make_combined_model_external_mets(org_list, model_name)

    lone_metabolites = [m for m in metabolite_exchange_model.metabolites if len(m.reactions) == 1]
    metabolite_exchange_model.remove_metabolites(lone_metabolites, destructive=True)

    return metabolite_exchange_model


# Function to save a json file in a d3flux readable format
def save_default_d3fluxmap_to_json(model_obj):

    #Create dictionary from cobra model
    default_map_dict = cobra.io.model_to_dict(model_obj)

    # Add in notes section dictionary
    # x,y coordinates will be added under map_info later
    for m in default_map_dict["metabolites"]:
        m['notes'] = {
            'map_info': {
                'display_name': m['id'][:-2].upper()
            }
        }
    # Save dict to json file
    default_map_json = json.dumps(default_map_dict, allow_nan=False)
    ex_met_exchange_map_fname = f"{model_obj.id}.json"
    with open(ex_met_exchange_map_fname, "w") as outfile:
        outfile.write(default_map_json)

    return default_map_dict


# Function to get external metabolites from a cobra model
# Note: assumes "_e" suffix represents an external metabolite
def get_ex_metabs(model_obj):
    org_m = [i.id for i in model_obj.metabolites]
    org_em = [i for i in org_m if i[-1] == 'e']

    return org_em


# Function to filter out metabolites by flux values
def filter_m_by_flux(org_ext_m, flux_df, exclude_0flux, flux_threshold):
    for i in org_ext_m:
        org_id = org_ext_m[i]['name']
        this_del_m_list = []

        # Loop through metabolites
        for m in org_ext_m[i]['ext_m']:
            keep_this_m = True
            this_m_name = "EX_" + m
            # Check if flux column exists for this metabolite
            if this_m_name in flux_df.columns:
                # Check for filtering criteria
                if np.isnan(flux_df.loc[org_id, this_m_name]):
                    keep_this_m = False
                elif exclude_0flux and (flux_df.loc[org_id, this_m_name] == 0):
                    keep_this_m = False
                elif (flux_threshold is not None) and (abs(flux_df.loc[org_id, this_m_name]) < flux_threshold):
                    keep_this_m = False

            # If filter criteria not met, add to removal list
            if not keep_this_m:
                this_del_m_list.append(m)

        # Update metabolite list for this organism
        org_ext_m[i]['ext_m'] = [m for m in org_ext_m[i]['ext_m'] if m not in this_del_m_list]

    return org_ext_m


# Function to re-organize external metabolites
def ex_metab_reorg(org_list, flux_df, exclude_0flux, flux_threshold):
    # Collect external metabolites from each of the three organisms
    org_ext_m = {
        1: {
            'name': org_list[0].id,
            'ext_m': get_ex_metabs(org_list[0]),
        },
        2: {
            'name': org_list[1].id,
            'ext_m': get_ex_metabs(org_list[1]),
        },
        3: {
            'name': org_list[2].id,
            'ext_m': get_ex_metabs(org_list[2]),
        }
    }

    # Remove metabolites for reactions according to flux map and user defined criteria
    if flux_df is not None:
        org_ext_m = filter_m_by_flux(org_ext_m, flux_df, exclude_0flux,
                                     flux_threshold)
    m1 = set(org_ext_m[1]['ext_m'])
    m2 = set(org_ext_m[2]['ext_m'])
    m3 = set(org_ext_m[3]['ext_m'])

    # Two-way shared external metabolites
    twoway_ext_m = {}
    twoway_ext_m["1_2"] = sorted(m1.intersection(m2).difference(m3))
    twoway_ext_m["1_3"] = sorted(m1.intersection(m3).difference(m2))
    twoway_ext_m["2_3"] = sorted(m2.intersection(m3).difference(m1))

    # Three-way shared external metabolites
    threeway_ext_m = sorted(m1.intersection(m2).intersection(m3))

    ### Create dict of shared external metabolites based on length and group type
    shared_ext_m = {'longest': {}, 'midlength': {}, 'shortest': {}}
    n_metabolites = {i: len(twoway_ext_m[i]) for i in twoway_ext_m.keys()}
    n_metabolites = list(
        reversed(sorted(n_metabolites.items(), key=lambda item: item[1])))

    def add_output(length, index_value):
        i, j = [int(n) for n in index_value.split('_')]
        shared_ext_m[length]['2way_key'] = index_value
        shared_ext_m[length]['orgs'] = [i, j]
        shared_ext_m[length]['org_names'] = [org_ext_m[i]['name'],
                                             org_ext_m[j]['name']]
        shared_ext_m[length]['ext_m'] = twoway_ext_m[index_value]

    add_output('longest', n_metabolites[0][0])
    add_output('midlength', n_metabolites[1][0])
    add_output('shortest', n_metabolites[-1][0])
    # add three-way list
    shared_ext_m['threeway'] = {
        'orgs': [1, 2, 3],
        'org_names': [org_ext_m[1]['name'], org_ext_m[2]['name'],
                      org_ext_m[3]['name']],
        'ext_m': threeway_ext_m
    }

    return shared_ext_m

# Function to assign positions for organism nodes
def get_org_pos(shared_ext_m, x_min, x_max, y_min, y_max):
    org_pos = {
        'longest_shortest': {  # bottom edge-left
            'name': [org for org in shared_ext_m['longest']['org_names'] if org in shared_ext_m['shortest']['org_names']][0],
            'x': x_min + (x_max - x_min) * (1.0 / 3.0),
            'y': y_max
        },
        'longest_midlength': {  # top edge-right
            'name': [org for org in shared_ext_m['longest']['org_names'] if org in shared_ext_m['midlength']['org_names']][0],
            'x': x_min + (x_max - x_min) * (2.0 / 3.0),
            'y': y_min
        },
        'shortest_midlength': {  # right edge-bottom
            'name': [org for org in shared_ext_m['midlength']['org_names'] if org in shared_ext_m['shortest']['org_names']][0],
            'x': x_max,
            'y': y_min + (y_max - y_min) * (2.0 / 3.0)
        }
    }
    return org_pos


def get_met_pos_limits(shared_ext_m, org_pos, x_min, x_max, y_min, y_max, d_x,
                       d_y):
    # Dict for the number of rows and columns for each list of metabolites
    row_cols = {}

    # Calculate about how many vertical space allocations go into a horizontal
    # one
    num_r2_cwidth = int(d_x/d_y)

    # The rest of this function will build a dict for the limits for each list
    # of metabolites .
    # Each group of values will need to be determined sequentially

    # Start with known limits
    pos_limits = {
        # top-left corner
        'longest': {
            # top-left corner
            'top_row_start': {
                'x': x_min,
                'y': y_min
            },
            # one column left of longest_midlength organism; one row above of
            # longest_shortest organism
            'bottom_row_end_max': {

                'x': org_pos['longest_midlength']['x'] - 0.5 * d_x,
                'y': org_pos['longest_shortest']['y'] - 0.5 * d_y
            }
        },
        # right edge - top
        'midlength': {
            # right side; one row below longest_midlength organism
            'top_row_end_min': {
                'x': x_max,
                'y': org_pos['longest_midlength']['y'] + d_y,
            },
            'bottom_row_end_max': {  # right side; one row above shortest_midlength organism
                'x': x_max,
                'y': org_pos['shortest_midlength']['y'] - d_y,
            },
            'right_col_center': {  # right edge-top
                'x': x_max,
                'y': y_min + (y_max - y_min) * (1.0 / 3.0)
            }
        },
        'shortest': {  # bottom edge - right
            # bottom left side; one column right of longest_shortest organism
            'bottom_row_start_min': {
                'x': org_pos['longest_shortest']['x'] + d_x,
                'y': y_max - d_y
            },
            # bottom right side; one column left of shortest_midlength organism
            'bottom_row_end_max': {
                'x': org_pos['shortest_midlength']['x'],
                'y': y_max
            },
            # bottom_edge right
            'bottom_row_center': {
                'x': x_min + (x_max - x_min) * (2.0 / 3.0),
                'y': y_max - d_y
            }
        }
    }

    # Calculate values for the Longest two-way list: top left corner
    L = len(shared_ext_m['longest']['ext_m'])

    row_col_max = {
        'longest': {
            'r_max': int(
                (pos_limits['longest']['bottom_row_end_max']['y'] - pos_limits['longest']['top_row_start']['y']) / d_y),
            'c_max': int(
                (pos_limits['longest']['bottom_row_end_max']['x'] - pos_limits['longest']['top_row_start']['x']) / d_x)
        }
    }
    if L > (row_col_max['longest']['r_max'] * row_col_max['longest']['c_max']):
        print("WARNING: The LONGEST list of metabolites is too long for the "
              "limits set in place")
    row_cols['longest'] = {}
    row_cols['longest']['minR'] = num_r2_cwidth * (1 + int(math.sqrt(L / num_r2_cwidth))) if num_r2_cwidth <= L else L
    row_cols['longest']['numC'] = max(int(L / row_cols['longest']['minR']), row_col_max['longest']['c_max'])#min(int(L / row_cols['longest']['minR']), row_col_max['longest']['c_max'])
    row_cols['longest']['numR'] = math.ceil(L / row_cols['longest']['numC'])

    ## Update limits based on longest list row/column info
    # Bottom right end of longest list placement position
    pos_limits['longest']['bottom_row_end'] = {
        'x': d_x * (row_cols['longest']['numC'] - 1) + pos_limits['longest']['top_row_start']['x'],
        'y': d_y * (row_cols['longest']['numR'] - 1) + pos_limits['longest']['top_row_start']['y']
    }
    # below longest list AND (right of longest list OR older x-limit)
    pos_limits['shortest']['top_row_start_min'] = {
        'x': max(pos_limits['shortest']['bottom_row_start_min']['x'],
                 d_x + pos_limits['longest']['bottom_row_end']['x']),
        'y': d_y + pos_limits['longest']['bottom_row_end']['y']
    }
    pos_limits['shortest']['bottom_row_start_min']['x'] = \
        pos_limits['shortest']['top_row_start_min']['x']

    pos_limits['threeway'] = {
        # right of longest list; below longest_midlength organism
        'top_row_start_min': {
            'x': d_x + pos_limits['longest']['bottom_row_end']['x'],
            'y': d_y + org_pos['longest_midlength']['y']
        }
    }

    # Calculate values for the Shortest two-way list:
    # centered on bottom-side right position
    L = len(shared_ext_m['shortest']['ext_m'])

    row_col_max['shortest'] = {
        'r_max': abs(int((pos_limits['shortest']['bottom_row_end_max']['y'] -
                          pos_limits['shortest']['top_row_start_min']['y']) / d_y)),
        'c_max': abs(int((pos_limits['shortest']['bottom_row_end_max']['x'] -
                          pos_limits['shortest']['top_row_start_min']['x']) / d_x))
    }

    if L > (row_col_max['shortest']['r_max'] * row_col_max['shortest']['c_max']):
        print("WARNING: The SHORTEST list of metabolites is too long for the "
              "limits set in place")

    row_cols['shortest'] = {}
    row_cols['shortest']['minC'] = (1 + int(math.sqrt(L / num_r2_cwidth))) if num_r2_cwidth <= L else int(math.sqrt(L))
    row_cols['shortest']['numR'] = min(int(L / row_cols['shortest']['minC']), row_col_max['shortest']['r_max'])
    row_cols['shortest']['numC'] = math.ceil(L / row_cols['shortest']['numR'])

    ## Update limits based on shortest list row/column info

    # Over-ride center position if list width goes beyond limit on either side
    radius = 0.5 * d_x * row_cols['shortest']['numC']
    max_right_width = abs(
        pos_limits['shortest']['bottom_row_end_max']['x'] - pos_limits['shortest']['bottom_row_center']['x'])
    max_left_width = abs(
        pos_limits['shortest']['top_row_start_min']['x'] - pos_limits['shortest']['bottom_row_center']['x'])

    if radius > max_right_width or radius > max_left_width:
        print("NOTE: Adjusting SHORTEST list center position away "
              "from originally defined location")
        pos_limits['shortest']['bottom_row_center']['x'] = \
            pos_limits['shortest']['top_row_start_min']['x'] + 0.5 * \
            abs(pos_limits['shortest']['bottom_row_end_max']['x'] -
                pos_limits['shortest']['top_row_start_min']['x'])

    # left of bottom-side right position
    pos_limits['shortest']['bottom_row_start'] = {
        'x': pos_limits['shortest']['bottom_row_center']['x'] - radius,
        'y': pos_limits['shortest']['bottom_row_center']['y']
    }

    # right of bottom-side right position
    pos_limits['shortest']['bottom_row_end'] = {
        'x': pos_limits['shortest']['bottom_row_center']['x'] + radius,
        'y': pos_limits['shortest']['bottom_row_center']['y']
    }
    if pos_limits['shortest']['bottom_row_start']['x'] < \
            pos_limits['shortest']['bottom_row_start_min']['x']:
        print("WARNING: The SHORTEST list of metabolites will overlap with "
              "the LONGEST list of metabolites")

    if pos_limits['shortest']['bottom_row_end']['x'] > \
            pos_limits['shortest']['bottom_row_end_max']['x']:
        print("WARNING: The SHORTEST list of metabolites goes beyond the "
              "horizontal limit position from the right")

    # directly above bottom_row_start position
    pos_limits['shortest']['top_row_start'] = {
        'x': pos_limits['shortest']['bottom_row_start']['x'],
        'y': pos_limits['shortest']['bottom_row_start']['y'] - d_y * (row_cols['shortest']['numR'] - 1)
    }
    if pos_limits['shortest']['top_row_start']['x'] < pos_limits['longest']['bottom_row_end']['x'] \
        or pos_limits['shortest']['top_row_start']['y'] < pos_limits['longest']['bottom_row_end']['y']:
        print("WARNING: The SHORTEST list of metabolites may overlap with "
              "the LONGEST list of metabolites")

    # directly above bottom_row_end position
    pos_limits['shortest']['top_row_end'] = {
        'x': pos_limits['shortest']['bottom_row_end']['x'],
        'y': pos_limits['shortest']['top_row_start']['y']
    }

    pos_limits['midlength']['bottom_row_end_max']['y'] = \
        min(
            pos_limits['midlength']['bottom_row_end_max']['y'],
            (pos_limits['shortest']['top_row_end']['y'] - d_y)
        )

    # Calculate values for the Midlength two-way list:
    # centered on right-side top position
    L = len(shared_ext_m['midlength']['ext_m'])

    row_col_max['midlength'] = {
        'r_max': abs(int((pos_limits['midlength']['bottom_row_end_max']['y'] -
                          pos_limits['midlength']['top_row_end_min']['y']) / d_y)),
        'c_max': 2
    }
    if L > (row_col_max['midlength']['r_max'] * row_col_max['midlength']['c_max']):
        print("WARNING: The MIDLENGTH list of metabolites is too long for the "
              "limits set in place")

    row_cols['midlength'] = {}
    row_cols['midlength']['numC'] = min(
        (1 + int(L / row_col_max['midlength']['r_max'])),
        row_col_max['midlength']['c_max']
    )
    row_cols['midlength']['numR'] = math.ceil(L / row_cols['midlength']['numC'])

    # Update limits based on midlength list row/column info

    # Over-ride center position if list width goes beyond limit on either side
    radius = 0.5 * d_y * row_cols['midlength']['numC']
    max_top_width = abs(
        pos_limits['midlength']['bottom_row_end_max']['y'] - pos_limits['midlength']['right_col_center']['y'])
    max_bottom_width = abs(
        pos_limits['midlength']['top_row_end_min']['y'] - pos_limits['midlength']['right_col_center']['y'])

    if radius > max_top_width or radius > max_bottom_width:
        print("NOTE: Adjusting MIDLENGTH list center position "
              "away from originally defined location")
        pos_limits['midlength']['right_col_center']['y'] = \
            pos_limits['midlength']['top_row_end_min']['y'] + \
            0.5 * abs(pos_limits['midlength']['bottom_row_end_max']['y']
                      - pos_limits['midlength']['top_row_end_min']['y'])

    pos_limits['midlength']['top_row_end'] = {
        'x': pos_limits['midlength']['top_row_end_min']['x'],
        'y': pos_limits['midlength']['right_col_center']['y'] - 0.5 * d_y * row_cols['midlength']['numR']
    }
    pos_limits['midlength']['bottom_row_end'] = {
        'x': pos_limits['midlength']['top_row_end']['x'],
        'y': pos_limits['midlength']['right_col_center']['y'] + 0.5 * d_y * row_cols['midlength']['numR']
    }
    pos_limits['midlength']['top_row_start'] = {
        'x': pos_limits['midlength']['top_row_end']['x'] - d_x * (row_cols['midlength']['numC'] - 1),
        'y': pos_limits['midlength']['top_row_end']['y']
    }
    pos_limits['midlength']['bottom_row_start'] = {
        'x': pos_limits['midlength']['top_row_start']['x'],
        'y': pos_limits['midlength']['bottom_row_end']['y']
    }

    pos_limits['threeway']['bottom_row_end_max'] = {
        'x': pos_limits['midlength']['top_row_start']['x'],
        'y': pos_limits['shortest']['bottom_row_start']['y']
    }

    ### Calculate values for the three-way list: center of figure
    L = len(shared_ext_m['threeway']['ext_m'])

    row_col_max['threeway'] = {
        'r_max': abs(int((pos_limits['threeway']['bottom_row_end_max']['y'] -
                          pos_limits['threeway']['top_row_start_min']['y']) / d_y)),
        'c_max': abs(int((pos_limits['threeway']['bottom_row_end_max']['x'] -
                          pos_limits['threeway']['top_row_start_min']['x']) / d_x))
    }

    if L > (row_col_max['threeway']['r_max'] * row_col_max['threeway']['c_max']):
        print("WARNING: The THREE-WAY list of metabolites is too long for the"
              " limits set in place")

    row_cols['threeway'] = {}
    row_cols['threeway']['minR'] = num_r2_cwidth * (1 + int(math.sqrt(L / num_r2_cwidth))) if num_r2_cwidth <= L else L
    row_cols['threeway']['numC'] = min(int(round(L / row_cols['threeway']['minR'], 0)),
                                       row_col_max['threeway']['c_max'])
    row_cols['threeway']['numR'] = math.ceil(L / row_cols['threeway']['numC'])

    # Identify center position for three-way list of metabolites
    pos_limits['threeway']['center'] = {
        'x': pos_limits['threeway']['top_row_start_min']['x'] +
             0.5 * abs(pos_limits['threeway']['bottom_row_end_max']['x'] -
                 pos_limits['threeway']['top_row_start_min']['x']),
        'y': pos_limits['threeway']['top_row_start_min']['y'] + 0.5 * abs(
            pos_limits['threeway']['bottom_row_end_max']['y'] - pos_limits['threeway']['top_row_start_min']['y'])
    }
    pos_limits['threeway']['top_row_start'] = {
        'x': pos_limits['threeway']['center']['x'] - 0.5 * d_x * row_cols['threeway']['numC'],
        'y': pos_limits['threeway']['center']['y'] - 0.5 * d_y * row_cols['threeway']['numR']
    }
    pos_limits['threeway']['top_row_end'] = {
        'x': pos_limits['threeway']['center']['x'] + 0.5 * d_x * row_cols['threeway']['numC'],
        'y': pos_limits['threeway']['center']['y'] - 0.5 * d_y * row_cols['threeway']['numR']
    }
    pos_limits['threeway']['bottow_row_start'] = {
        'x': pos_limits['threeway']['center']['x'] - 0.5 * d_x * row_cols['threeway']['numC'],
        'y': pos_limits['threeway']['center']['y'] + 0.5 * d_y * row_cols['threeway']['numR']
    }
    pos_limits['threeway']['bottow_row_end'] = {
        'x': pos_limits['threeway']['center']['x'] + 0.5 * d_x * row_cols['threeway']['numC'],
        'y': pos_limits['threeway']['center']['y'] + 0.5 * d_y * row_cols['threeway']['numR']
    }

    return row_cols, pos_limits


# Function to create dataframe of node coordinate positions
def get_node_pos(org_pos, shared_ext_m, pos_limits, row_cols, d_x, d_y):

    # Create empty df with column names
    ext_met_df = pd.DataFrame(columns=['Type', 'DisplayName', 'x', 'y'])

    # Include organism nodes into df
    for i in org_pos:
        ext_met_df.loc[len(ext_met_df.index)] = ['organism', f"{org_pos[i]['name']}_", org_pos[i]['x'], org_pos[i]['y']]

    # Include metabolite nodes into df
    for k in shared_ext_m.keys():
        this_start_x = pos_limits[k]['top_row_start']['x']
        this_start_y = pos_limits[k]['top_row_start']['y']
        this_r = row_cols[k]['numR']
        this_c = row_cols[k]['numC']
        row_count = 0
        col_count = 0
        for m in shared_ext_m[k]['ext_m']:
            this_x = this_start_x + col_count * d_x
            this_y = this_start_y + row_count * d_y
            ext_met_df.loc[len(ext_met_df.index)] = ['metabolite', m[:-2].upper(), this_x, this_y]
            row_count += 1
            if row_count == this_r:
                row_count = 0
                col_count += 1
            if row_count == this_r and col_count == this_c:
                print("ERROR: went over column limit")

    # Set displayname as df index
    ext_met_df = ext_met_df.set_index('DisplayName')

    return ext_met_df


# Function to translate the flux_df into a useful flux_dict
def create_flux_dict_info(model_dict, flux_df, exclude_0flux=False,
                          flux_threshold=0.0):

    # flattens all columns. Can now filter it all easily
    melted = pd.melt(flux_df.reset_index(), id_vars=['index'])
    melted['rxn_ids'] = melted['variable'] + '_' + melted['index']
    # filter to only reactions in model
    rxns = [r['id'] for r in model_dict['reactions']]
    valid_rxns = melted.loc[melted.rxn_ids.isin(rxns)].copy()

    # only columns we need
    cols = ['variable', 'rxn_ids', 'value']
    # remove nans
    rxns_to_remove = set(valid_rxns[valid_rxns.value.isna()]['rxn_ids'].values)

    flux = valid_rxns[~valid_rxns.isna()][cols].set_index('variable')

    if flux_threshold or exclude_0flux:
        too_small = flux.loc[np.abs(flux['value']) <= flux_threshold]
        flux = flux.loc[np.abs(flux['value']) > flux_threshold]
        # remove values too small
        rxns_to_remove = rxns_to_remove.union(set(too_small.rxn_ids))

    excluded_metabolites = set(rxn[3:-6] for rxn in rxns_to_remove)

    flux_values = flux.to_dict()['value']
    flux_dict_info = {
        'flux_dict': flux_values,
        'excluded_met': list(excluded_metabolites),
        'excluded_reacts': list(rxns_to_remove)
    }

    return flux_dict_info


# Function to generate d3flux flux map
def create_flux_map(model_name, ext_met_df, flux_dict=None):
    # Read in JSON file using expected file naming convention
    ex_met_exchange_map_fname = f"{model_name}.json"
    with open(ex_met_exchange_map_fname, "r") as read_file:
        ex_met_exchange_map = json.load(read_file)

    # Update x-y coordinates for metabolites in the JSON object
    for m in ex_met_exchange_map['metabolites']:

        # m_name = m['notes']['map_info']['display_name']
        # Check if m is for one of the three organisms for appropriate m_name choice
        m_name = f"{m['id']}_" if m['id'] == m['compartment'] \
            else m['notes']['map_info']['display_name']

        if m_name in ext_met_df.index:
            m['notes']['map_info']['x'] = ext_met_df.loc[m_name]['x']
            m['notes']['map_info']['y'] = ext_met_df.loc[m_name]['y']
    json_object = json.dumps(ex_met_exchange_map)

    # Create updated filename to prevent override
    updated_fname = f"{model_name}_Updated.json"
    f_count = 0
    while os.path.exists(updated_fname):
        f_count += 1
        updated_fname = f"{model_name}_Updated_{f_count}.json"

    # Print to screen the filename used to generate this figure
    print("Showing map corresponding to the following file: " + updated_fname)

    # Write JSON object to file
    with open(updated_fname, "w") as outfile:
        outfile.write(json_object)

    # Produce flux map model and object
    fluxmap_dict = {}
    fluxmap_dict['fname'] = updated_fname
    fluxmap_dict['model'] = cobra.io.load_json_model(fluxmap_dict['fname'])
    fluxmap_dict['map'] = d3f.flux_map(fluxmap_dict['model'], flux_dict=flux_dict)

    # fluxmap_dict_model = cobra.io.load_json_model(updated_fname)
    # fluxmap_dict_map = d3f.flux_map(fluxmap_dict_model)

    # ex_met_exchange_model_new = cobra.io.load_json_model(updated_fname)
    # fluxmap_viz = d3f.flux_map(ex_met_exchange_model_new)

    # return fluxmap_dict_model, fluxmap_dict_map
    return fluxmap_dict


# Function to produce reorganized visualization for shared external metabolites
# for 3-organisms
def make_3organism_extmetab_viz(org_list, model_name,
                                multi_org_model=None, reorg_after_filter=False,
                                flux_df=None, exclude_0flux=False, flux_threshold=None,
                                x_min=0, x_max=800, d_x=120,
                                y_min=0, y_max=650, d_y=20):
    # INPUTS:
    # org_list = a list of three cobra models, one for each organism
    # model_name = model name used for combined model id name and corresponding json filenames
    # multi_org_model =
    # flux_df =
    # exclude_0flux
    # flux_threshold
    # x_min, x_max, y_min, y_max = min & max limits for x & y coordinates for placing nodes in final figure
    # d_x, d_y = minimum spacing between nodes in the x & directions

    # Create a multi-organism model if one isn't already provided
    if multi_org_model is None:
        multi_org_model = make_combined_model_external_mets_shared_only(
            org_list, model_name
        )

    # Save multi-organism model in JSON format to work with
    model_dict = save_default_d3fluxmap_to_json(multi_org_model)

    # Reorganize external metabolites shared between 2 and 3 organisms
    # TODO: output list of excluded metabolites corresponding to flux criteria
    # Any way to cross reference this with the multi-org model?
    # Careful not to omit the original organisms which are represented
    # as metabolite nodes

    # There is no difference between function calls or args - JP
    if reorg_after_filter:
        # Filter out using flux data first, then reorganize the node coordinates
        shared_ext_m = ex_metab_reorg(
            org_list, flux_df, exclude_0flux, flux_threshold
        )
    else:
        # First create model/dict using all metabolites, then generate the
        # d3flux map using exclusion info (later)
        shared_ext_m = ex_metab_reorg(
            org_list, flux_df, exclude_0flux, flux_threshold
        )

    # Assign positions for single nodes for the three organisms
    org_pos = get_org_pos(shared_ext_m, x_min, x_max, y_min, y_max)

    # Get position limits for each metabolite node
    row_cols, pos_limits = get_met_pos_limits(
        shared_ext_m, org_pos, x_min, x_max, y_min, y_max, d_x, d_y
    )

    # Get coordinate positions for each external metabolite node,
    # including organism nodes
    ext_met_df = get_node_pos(
        org_pos, shared_ext_m, pos_limits, row_cols, d_x, d_y
    )

    # There is no difference between function calls or args - JP
    # Create the flux map

    fluxmap_dict = create_flux_map(model_name, ext_met_df)

    return fluxmap_dict


    # if reorg_after_filter:
    #     fluxmap_dict = create_flux_map(model_name, ext_met_df)
    #     flux_dict_info = create_flux_dict_info(model_dict, flux_df,
    #                                            exclude_0flux, flux_threshold)
    # else:
    #     fluxmap_dict = create_flux_map(model_name, ext_met_df)
    #     flux_dict_info = create_flux_dict_info(model_dict, flux_df,
    #                                            exclude_0flux, flux_threshold)
    #
    # return fluxmap_dict, flux_dict_info


# Save d3flux readable model to json file
# save_default_d3fluxmap_to_json(av_syn_rt_exchange_model)


