import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objects as go


def create_fluxflow_viz(consort_model, flux_df, org2group_dict,
                        flux_threshold=0.01):

    x_pos = {
        "source_m": 0.001,  # source metabolites that go into the system
        "donor_org": 0.25,  # organisms acting as donors
        "shared_m": 0.5,  # metabolites shared between more than two organisms
        "receiver_org": 0.75,  # organisms acting as receivers
        "excreted_m": 0.999,  # output metabolites that go into the system
    }
    # Update flux_df by filtering to external metabs and refining column names
    f = refine_flux_df(flux_df)

    # Get metabolite carrier groupings, including user definitions for organism
    carrier_m, carrier_names = get_metab_carrier_groupings(consort_model, org2group_dict)

    # Get dataframe with node information
    node_df, links_df = get_nodes_and_links(x_pos, f, carrier_m, carrier_names, flux_threshold)

    # Generate diagram object
    viz_obj = get_diagram(node_df, links_df)

    return viz_obj


def refine_flux_df(flux_df):
    org_names = get_org_names(flux_df)
    f = flux_df.replace(0, np.nan)
    good = (f.isnull().sum()<len(org_names)+1)
    good = list(good[good].index.values)
    f = f[good].copy()
    # Get list of external metabolites
    ex_metabs_cols = [c for c in f.columns if
                      (c.startswith('EX_') and c.endswith('_e'))]
    ex_metabs_names = [m[3:-2] for m in ex_metabs_cols]

    # Create df for only external metabolites (renamed) and organisms
    ex_flux_df = f.loc[org_names][ex_metabs_cols]
    ex_flux_df.columns = ex_metabs_names

    return ex_flux_df


def get_metab_carrier_groupings(model, org2group_dict):
    # Create full metabolite carrier groupings
    carrier_m = initialize_carrier_m_groupings()

    # Add organisms to user selected mtabolite group
    for k, v in org2group_dict.items():
        carrier_m[v]['names'].append(k)

    # Add metabolites corresponding to their elemental composition
    for m in model.metabolites:
        if m.id.lower().startswith("photon"):
            carrier_m['sunlight']['names'].append(m.id[:-2])
        elif ('C' in m.elements) and ('N' not in m.elements):
            carrier_m['C']['names'].append(m.id[:-2])
        elif ('C' not in m.elements) and ('N' in m.elements):
            carrier_m['N']['names'].append(m.id[:-2])
        elif ('C' in m.elements) and ('N' in m.elements):
            carrier_m['C+N']['names'].append(m.id[:-2])
        elif m.id in carrier_m['jet_fuel']['list_of_m']:
            carrier_m['jet_fuel']['names'].append(m.id[:-2])

    # Get a separate variable for metabolites that have been grouped
    carrier_names = []
    for k in carrier_m:
        carrier_names.extend(carrier_m[k]['names'])

    return carrier_m, carrier_names


def get_org_names(flux_df):
    # Get list of org names from flux_df
    return sorted([org for org in flux_df.index if not(org=='medium')])


def get_ext_metab_names(flux_df):
    return [c[3:-2] for c in flux_df.columns if (c.startswith("EX_") and c.endswith("_e"))]


def get_nodes_and_links(x_pos, flux_df, carrier_m, carrier_names,
                        flux_threshold=0.01):
    # Initialize node df
    node_df = initialize_node_df(flux_df, x_pos)

    # Refine node_df with links between nodes and update x_pos's
    node_df, m_in_use, links_df = refine_node_df(node_df, flux_df, x_pos,
                                                 flux_threshold)

    # Assign colors to nodes
    org_names = get_org_names(flux_df)
    node_df = assign_colors_node_df(node_df, org_names, m_in_use, carrier_m,
                                    carrier_names)

    # Assign y-positions to nodes based on color & x-position groupings
    node_df = get_y_pos(node_df, carrier_m)

    return node_df, links_df


def initialize_node_df(flux_df, x_pos):

    # Get organism name info
    org_names = get_org_names(flux_df)
    n_orgs = len(org_names)

    # Add Donor organisms to node df
    doner_df = pd.DataFrame.from_dict({
        "Node_Name": org_names,
        "Type": ["donor_org"]*n_orgs,
        "Xpos": [x_pos['donor_org']]*n_orgs,
        "Ypos": [np.nan]*n_orgs,
        "Color": ['rgb(0,0,0)']*n_orgs
    })

    # Add Receiver organisms to node df
    receiver_df = pd.DataFrame.from_dict({
        "Node_Name": org_names,
        "Type": ["receiver_org"]*n_orgs,
        "Xpos": [x_pos['receiver_org']]*n_orgs,
        "Ypos": [np.nan]*n_orgs,
        "Color": ['rgb(0,0,0)']*n_orgs
    })

    # Get external metabolite names from flux_df
    ex_metabs_names = flux_df.columns

    # Add external metabolites to node df
    # Note: needs refinement after considering reactions (see refine_node_df())
    Num_ms = len(ex_metabs_names)

    external_df = pd.DataFrame.from_dict({
        "Node_Name": ex_metabs_names,
        "Type": ["metabolites"]*Num_ms, # Will need updating
        "Xpos": [np.nan]*Num_ms, # Will need updating
        "Ypos": [np.nan]*Num_ms, # Will need updating,
        "Color": ['rgb(0,0,0)']*Num_ms
    })
    node_df = pd.concat([doner_df, receiver_df, external_df], ignore_index=True)
    node_df.reset_index(names=['Node_Idx'], inplace=True)
    return node_df


def refine_node_df(node_df, flux_df, x_pos, flux_threshold):
    # Get organism name info
    org_names = get_org_names(flux_df)
    n_orgs = len(org_names)

    # Generate links between nodes and refine x-positions of nodes
    # Construct link df - first pass at connecting organisms to metabolites
    output = []

    m_in_use = set()
    name_to_index = node_df.set_index('Node_Name').to_dict()['Node_Idx']

    name_to_type_index = node_df.groupby('Node_Name')[['Type', 'Node_Idx']].apply(
        lambda x: x.set_index('Type').to_dict(orient='index')
    ).to_dict()

    for this_m in flux_df.columns:

        this_met_idx = name_to_index[this_m]
        this_col = flux_df.loc[org_names, this_m]
        # Count good fluxes and flow direction
        sum_neg_flux = ((np.abs(this_col) > flux_threshold) & (this_col < 0)).sum()
        sum_pos_flux = ((np.abs(this_col) > flux_threshold) & (this_col > 0)).sum()
        sum_good_flux = n_orgs - ((np.abs(this_col) < flux_threshold) | (np.isnan(this_col))).sum()
        if sum_good_flux < 1:
            continue
        # Include metabolite into tracked list
        m_in_use.add(this_m)

        if sum_pos_flux == 0:
            node_type =  'source_m'
            x_position = x_pos['source_m']
        elif sum_neg_flux == 0:
            node_type =  'excreted_m'
            x_position = x_pos['excreted_m']
        else:
            node_type =  'shared_m'
            x_position = x_pos['shared_m']
            spacing = 0.15
            # More excreted than sourced
            if sum_pos_flux > sum_neg_flux:
                x_position -= spacing
            # More sourced than excreted
            else:
                x_position += spacing

        node_df.loc[this_met_idx, "Type"] = node_type
        node_df.loc[this_met_idx, "Xpos"] = x_position

        for this_org in org_names:
            this_val = flux_df.loc[this_org, this_m]
            # below threshold
            if np.isnan(this_val) or abs(this_val) < flux_threshold:
                continue

            # Reaction from media to donor-org: only - flux
            if sum_pos_flux == 0:
                this_org_idx = name_to_type_index[this_org]['donor_org']['Node_Idx']
                output.append([this_met_idx, this_m, this_org_idx, this_org, abs(this_val)])

            # Reaction from receiver-org to excreted-met: overall, only + flux
            elif sum_neg_flux==0:
                this_org_idx = name_to_type_index[this_org]['receiver_org']['Node_Idx']
                output.append([ this_org_idx, this_org, this_met_idx, this_m, abs(this_val)])

            # Reaction between shared-met and an org
            elif sum_pos_flux > 0 and sum_neg_flux > 0:
                # + means secretion from the organism into the environment
                if this_val > 0:
                    this_org_idx = name_to_type_index[this_org]['donor_org']['Node_Idx']
                    output.append([this_org_idx, this_org, this_met_idx, this_m, abs(this_val)])

                # - means uptake into the organism from the environment
                elif this_val < 0:
                    this_org_idx = name_to_type_index[this_org]['receiver_org']['Node_Idx']
                    output.append([this_met_idx, this_m, this_org_idx, this_org, abs(this_val)])

            else:
                print("No option from above")

    links_df = pd.DataFrame(output, columns=["Sources", "Source_Name", "Targets", "Target_Name", "Values"])
    return node_df, m_in_use, links_df


def get_y_pos(node_df, carrier_m):
    # Carrier type ordering
    c_type_order = ['sunlight', 'C', 'N', 'C+N', 'jet_fuel']
    for i, group in node_df.groupby('Type'):

        this_n_list = list(group['Node_Name'])
        this_yPos_N = len(this_n_list)
        # Add a spacer at the start of each group
        this_yPos_N += 1
        new_n_list = ['space']

        # Collect node-names in preferred ordered list
        for this_c_type in c_type_order:
            names_in_this_c_type = [n for n in group['Node_Name'] if
                                    n in carrier_m[this_c_type]['names']]

            if len(names_in_this_c_type) > 0:
                # Include these names in list
                new_n_list.extend(names_in_this_c_type)

                # Add a spacer after each group
                this_yPos_N += 1
                new_n_list.append('space')

                # Remove included nodes from this_n_list
                for n in names_in_this_c_type:
                    j = [i for i, this_n in enumerate(this_n_list) if this_n == n][0]
                    this_n_list.pop(j)

        # Collect remaining node-names into new_n_list
        new_n_list.extend(this_n_list)

        # Add ending space if it's missing
        if new_n_list[-1] != 'space':
            this_yPos_N += 1
            new_n_list.append('space')

        new_n_df = pd.DataFrame.from_dict(
            {'Name': new_n_list, 'Ypos': np.linspace(0.01, 0.9, this_yPos_N)})
        # Update node_df with y-values
        for i, n in enumerate(new_n_df['Name']):
            if not (n == 'space'):
                node_df.loc[node_df['Node_Name'] == n, 'Ypos'] = new_n_df.loc[i, 'Ypos']

    return node_df


def initialize_carrier_m_groupings():
    carrier_m = {
        'C': {
            'names': [],
            'color': 'rgb(78, 173, 91)' # Hex = 4EAD5B #'green'
        },
        'N': {
            'names': [],
            'color': 'rgb(202, 123, 44)' # Hex = CA7B2C # 'orange', #
        },
        'C+N': {
            'names': [],
            'color': 'rgb(0,0,1)' # 'blue',
        },
        'sunlight': {
            'names': [],
            'color': 'rgb(245, 194, 66)' # Hex = F5C242 # 'yellow',#
        },
        'jet_fuel': {
            'list_of_m': [], # NEEDS TO BE FILLED WITH JET-FUEL CANDIDATE METABOLITES
            'names': [],
            'color': 'rgb(133, 26, 43)' # Hex = 851A2B # 'red', #
        },
    }
    return carrier_m


def assign_colors_node_df(node_df, org_names, m_in_use, carrier_m,
                          carrier_names):
    # Color palette in hexadecimal format
    palette = np.array(px.colors.qualitative.Light24)

    # Color palette as RGB tuples
    palette_rgb = [tuple(int(c.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4))
                   for c in palette]
    palette_rgb_list = [f'rgb{c}' for c in palette_rgb]

    # Add colors to node_df
    color_count = 0
    org_idx = [i for i, n in enumerate(node_df['Node_Name']) if
               (n in org_names)]
    m_in_use_idx = [i for i, n in enumerate(node_df['Node_Name']) if
                    (n in m_in_use)]

    for i in org_idx + m_in_use_idx:
        this_node = node_df.loc[i, 'Node_Name']
        if this_node in carrier_names:
            for this_type in carrier_m:
                if this_node in carrier_m[this_type]['names']:
                    node_df.loc[i, 'Color'] = carrier_m[this_type]['color']
        else:
            node_df.loc[i, 'Color'] = palette_rgb_list[color_count]
            color_count += 1
            if color_count == len(palette_rgb_list):
                color_count = 0

    return node_df


def get_diagram(node_df, links_df):

    fig_obj = go.Figure(go.Sankey(
        valueformat=".3f",
        node=dict(
            pad=15,
            thickness=30,
            label=node_df['Node_Name'].values,
            x=node_df['Xpos'].values,
            # y=list(node_df['Ypos'].values),
            color=node_df['Color'].values
        ),
        link=dict(
            arrowlen=30,
            source=links_df["Sources"].values,
            target=links_df["Targets"].values,
            value=links_df["Values"].values
        )
    ))

    fig_obj.update_layout(width=1200, height=700, margin={'t': 10, 'b': 10})
    fig_obj.update_layout(font_size=16)

    return fig_obj