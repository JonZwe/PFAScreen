from collections import Counter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import networkx as nx
import pandas as pd
import pyopenms as oms
from matchms import Spectrum, calculate_scores
from matchms.similarity import CosineGreedy
from matchms.networking import SimilarityNetwork



def gnps_layout(G, 
                layout_func=nx.kamada_kawai_layout, 
                padding=0.3):
    
    """Grid-pack connected components (largest first), each laid out independently."""
    components = sorted(nx.connected_components(G), key=len, reverse=True)

    sub_layouts = []
    for comp in components:
        sub = G.subgraph(comp)
        if len(sub) == 1:
            pos = {next(iter(sub.nodes())): np.array([0.0, 0.0])}
        else:
            pos = layout_func(sub)

        coords = np.array(list(pos.values()))
        coords -= coords.mean(axis=0)
        span = coords.max(axis=0) - coords.min(axis=0)
        span[span == 0] = 1.0
        coords /= span

        scale = np.sqrt(len(sub))
        coords *= scale
        sub_layouts.append((dict(zip(pos.keys(), coords)), scale))

    pos_final = {}
    x_cursor, y_cursor, row_height = 0.0, 0.0, 0.0
    row_width = np.sqrt(sum(s for _, s in sub_layouts)) * 4

    for coords, scale in sub_layouts:
        box_size = scale * 2 + padding
        if x_cursor + box_size > row_width and x_cursor > 0:
            x_cursor, y_cursor, row_height = 0.0, y_cursor + row_height + padding, 0.0

        for node, xy in coords.items():
            pos_final[node] = xy + np.array([x_cursor + scale, y_cursor + scale])

        x_cursor += box_size
        row_height = max(row_height, box_size)

    return pos_final


def molecular_network(df, 
                      identifier_key='precursor_mz', 
                      col_highlight = 'compound_names',
                      colormap_column='rt',
                      score_cutoff=0.7, 
                      max_links=10) -> None:
    """
    Create a molecular network based in functionality of matchms
    Inspiration from: 
    https://blog.esciencecenter.nl/build-a-mass-spectrometry-analysis-pipeline-in-python-using-matchms-part-iii-molecular-91891248ee34
    """
    # get only all features with MS2 spectra
    df_ms2 = df[df['mzs_ms2'].notna()].copy()
    df_ms2['highlight'] = df_ms2[col_highlight].notna()
    highlighted = df_ms2[df_ms2['highlight'] == True]['mz'].round(4).tolist()

    # NOTE: Add all possible metadata!!
    spectra = []
    for n in range(len(df_ms2)):
        spectra.append(Spectrum(mz=np.array(df_ms2['mzs_ms2'].iloc[n]),
                                intensities=np.array(df_ms2['ints_ms2'].iloc[n]),    
                                metadata={"precursor_mz": df_ms2['mz'].iloc[n], 
                                          "precursor_rt": df_ms2['rt'].iloc[n],
                                          "spectrum_id": df_ms2.index[n]}))

    # calculate scores between all spectra to create a similarity matrix
    scores = calculate_scores(spectra, spectra, CosineGreedy())
    ms_network = SimilarityNetwork(identifier_key=identifier_key, score_cutoff=score_cutoff, max_links=max_links)
    ms_network.create_network(scores, 'CosineGreedy_score')
    my_network = ms_network.graph

    # Add rounded m/z as node attribute
    mz_to_rt = df_ms2.set_index('mz')['rt'].to_dict()

    for node in my_network.nodes:
        mz = float(node)
        my_network.nodes[node]['rounded_mz'] = round(mz, 4)
        my_network.nodes[node]['rt'] = mz_to_rt.get(mz, np.nan)

    fig, ax = plt.subplots(figsize=(10, 10))
    labels = {node: my_network.nodes[node]['rounded_mz'] for node in my_network.nodes}
    #pos = nx.spring_layout(my_network, seed=42)
    pos = gnps_layout(my_network)

    # Fixed fill color for all nodes
    if colormap_column == None:
        node_fill_color = 'darkcyan'
    else:
        node_fill_color = [my_network.nodes[node][colormap_column] for node in my_network.nodes]

    node_size = 250

    # Determine edge colors based on highlighting
    node_border_colors = []
    for node in my_network.nodes:
        if my_network.nodes[node]['rounded_mz'] in highlighted:
            node_border_colors.append('fuchsia')
        else:
            node_border_colors.append('grey')  # Or any other color you want for default borders

    # Draw nodes with different edge colors (node_color is the fill, edgecolors is the border)
    nx.draw_networkx_nodes(my_network, pos, node_size=node_size, node_color=node_fill_color, cmap=plt.cm.viridis, edgecolors=node_border_colors, linewidths=1.5)
    nx.draw_networkx_edges(my_network, pos, alpha=0.3, edge_color='blue')
    nx.draw_networkx_labels(my_network, pos, labels=labels, font_size=10, font_color='black', alpha=0.7)

    plt.title(f"Network: {len(df_ms2)} MS2 spectra | {len(my_network.nodes)} in network (with cosine > {score_cutoff})")
    plt.axis("off")
    if colormap_column == None:
        plt.show()
    else:
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(node_fill_color), vmax=max(node_fill_color)))
        fig.colorbar(sm, ax=ax, label=colormap_column)
        plt.show()
    
    #subnetworks = [my_network.subgraph(c).copy() for c in nx.connected_components(my_network)]
    #clusters = {}
    #for i, subnet in enumerate(subnetworks):
    #    if len(subnet.nodes) > 2: # exclude clusters with <=2 nodes
    #        clusters[i] = len(subnet.nodes)
    #net.plot_cluster(subnetworks[19])
    #nx.write_graphml(my_network, "network.graphml")



def single_molecular_network(df, 
                             feature_id,
                             identifier_key='precursor_mz', 
                             col_highlight = 'compound_names',
                             colormap_column='rt',
                             score_cutoff=0.7) -> None:
    # NOTE: DEV STAGE! DOES NOT WORK!

    """
    Create network for a single feature based in functionality of matchms
    """

    feature = df.loc[feature_id]
    spectrum_ref = Spectrum(mz=np.array(feature['mzs_ms2']), 
                            intensities=np.array(feature['ints_ms2']),   
                            metadata={"precursor_mz": feature['mz'], 
                                      "precursor_rt": feature['rt'],
                                      "spectrum_id": feature_id})

    # get only all features with MS2 spectra
    df_ms2 = df[df['mzs_ms2'].notna()].copy()
    df_ms2['highlight'] = df_ms2[col_highlight].notna()

    # NOTE: Add all possible metadata!!
    spectra = []
    for n in range(len(df_ms2)):
        spectra.append(Spectrum(mz=np.array(df_ms2['mzs_ms2'].iloc[n]),
                                intensities=np.array(df_ms2['ints_ms2'].iloc[n]),    
                                metadata={"precursor_mz": df_ms2['mz'].iloc[n], 
                                          "precursor_rt": df_ms2['rt'].iloc[n],
                                          "spectrum_id": df_ms2.index[n]}))

    # calculate scores between all spectra to create a similarity matrix
    scores = calculate_scores([spectrum_ref], spectra, CosineGreedy())
    cosine_scores = scores.to_array("CosineGreedy_score")[0]  # shape (N,)

    G = nx.Graph()

    # Add reference node with attributes
    G.add_node(feature_id,
            mz=feature['mz'],
            rounded_mz=round(feature['mz'], 4),
            rt=feature['rt'])

    # Add edges to all spectra above threshold
    for spec, score in zip(spectra, cosine_scores):
        if score > score_cutoff:
            node_id = spec.metadata["spectrum_id"]  # use spectrum_id as node ID
            G.add_node(node_id,
                    mz=spec.metadata["precursor_mz"],
                    rounded_mz=round(spec.metadata["precursor_mz"], 4),
                    rt=spec.metadata["precursor_rt"])
            G.add_edge(feature_id, node_id, weight=score)

    # -----------------------------
    # 4. Plot the network
    # -----------------------------
    plt.figure()

    #pos = nx.spring_layout(G, seed=42)
    pos = gnps_layout(G)
    labels = {node: G.nodes[node]["rounded_mz"] for node in G.nodes}

    nx.draw_networkx_nodes(
        G, pos,
        node_size=250,
        node_color='red',
        cmap=plt.cm.viridis,
        edgecolors='black',
        linewidths=1.5
    )

    # Draw edges
    nx.draw_networkx_edges(
        G, pos,
        alpha=0.3,
        edge_color="blue"
    )

    # Labels
    nx.draw_networkx_labels(
        G, pos,
        labels=labels,
        font_size=10,
        font_color="black",
        alpha=0.7
    )

    # Title and colorbar
    plt.title(f"Network: {len(G.nodes)} nodes | {len(G.edges)} edges (cosine > {score_cutoff})")
    plt.axis("off")

    plt.show()


def mass_difference_network(mz_array,
                            diffs=None,
                            mz_tol=0.005,
                            node_color=None,
                            node_color_label='m/z',
                            node_labels=None,
                            mz_col='mz',
                            node_size=300, 
                            font_size=6,
                            cmap='coolwarm',
                            show_group_number=False
                            ) -> tuple:
    """
    Create annotated mass difference network plots from m/z values and given mass differences.

    Parameters
    ----------
    mz_array : array-like or pandas.DataFrame
        m/z values, or a DataFrame containing an m/z column.
    diffs : list
        Mass differences as formulas or numeric values.
    mz_tol : float
        Mass tolerance for matching differences.
    network_n : int
        Rank of connected component to plot separately, sorted by component size.
    k : float
        Spring layout spacing parameter.
    node_color : array-like, pandas.Series, or str, optional
        Values used for node coloring. If mz_array is a DataFrame, this can be a column name.
        If omitted, nodes are colored by m/z.
    node_color_label : str, optional
        Label for the node colorbar.
    node_labels : array-like, pandas.Series, or str, optional
        Labels shown on nodes. If mz_array is a DataFrame, this can be a column name.
        If omitted, rounded m/z values are used.
    mz_col : str
        Column containing m/z values when mz_array is a DataFrame.
    """
    if diffs is None:
        diffs = ['CF2']
    cmap_nodes = plt.get_cmap(cmap)

    input_is_df = isinstance(mz_array, pd.DataFrame)

    if input_is_df:
        source_df = mz_array
        mz_values = source_df[mz_col].to_numpy()
    else:
        source_df = None
        mz_values = np.asarray(mz_array)

    def _resolve_node_values(values, value_name):
        if values is None:
            return None

        if input_is_df and isinstance(values, str):
            resolved = source_df[values].to_numpy()
        else:
            resolved = np.asarray(values)

        if len(resolved) != len(mz_values):
            raise ValueError(
                f"{value_name} must have the same length as mz_array "
                f"({len(resolved)} != {len(mz_values)})."
            )

        return resolved

    node_color_values = _resolve_node_values(node_color, 'node_color')
    node_label_values = _resolve_node_values(node_labels, 'node_labels')

    if node_color_values is None:
        node_color_values = mz_values
        if node_color_label is None:
            node_color_label = 'm/z'

    rounded_mz = np.round(mz_values, 4)

    # If rounded m/z values are duplicated, the graph nodes are merged by NetworkX.
    # The first matching row is used for node colors and labels.
    mz_to_idx = {}
    for idx, mz in enumerate(rounded_mz):
        mz_to_idx.setdefault(mz, idx)

    # Compute exact mass differences
    m_diffs = np.array([
        oms.EmpiricalFormula(diff).getMonoWeight() if isinstance(diff, str) else diff
        for diff in diffs
    ])

    # Compute difference matrix
    mzx, mzy = np.meshgrid(mz_values, mz_values)
    diff_matrix = abs(mzx - mzy)

    # Find matching indices
    connect_idx_1 = []
    connect_idx_2 = []
    diff_connection = []

    for n, m_diff in enumerate(m_diffs):
        log_matrix = np.logical_and(
            diff_matrix >= m_diff - mz_tol,
            diff_matrix <= m_diff + mz_tol
        )
        idx1, idx2 = np.where(log_matrix)

        # Keep only one direction for each pair.
        keep = idx1 < idx2
        idx1 = idx1[keep]
        idx2 = idx2[keep]

        connect_idx_1.extend(idx1)
        connect_idx_2.extend(idx2)
        diff_connection.extend([diffs[n]] * len(idx1))

    # Build DataFrame
    df = pd.DataFrame({
        'mass_1': rounded_mz[np.array(connect_idx_1)],
        'mass_2': rounded_mz[np.array(connect_idx_2)],
        'connect': diff_connection
    })

    diff_to_color = {diff: i for i, diff in enumerate(diffs)}
    df['color'] = df['connect'].map(diff_to_color)

    # Build NetworkX graph
    G = nx.from_pandas_edgelist(
        df,
        source='mass_1',
        target='mass_2',
        edge_attr=['color', 'connect']
    )

    if len(G.nodes) == 0:
        print("No mass-difference connections found.")
        return df, diff_matrix

    # Layout
    #pos = nx.spring_layout(G, k=k)
    pos = gnps_layout(G)
    for node, (x, y) in pos.items():
        G.nodes[node]['x'] = x
        G.nodes[node]['y'] = y

    # Node colors and labels
    graph_node_colors = [
        node_color_values[mz_to_idx[node]]
        for node in G.nodes()
    ]

    if node_label_values is None:
        graph_node_labels = {
            node: str(node)
            for node in G.nodes()
        }
    else:
        graph_node_labels = {
            node: str(node_label_values[mz_to_idx[node]])
            if pd.notna(node_label_values[mz_to_idx[node]]) else ''
            for node in G.nodes()
        }

    edge_colors = [a['color'] for _, _, a in G.edges(data=True)]
    edge_labels = {(u, v): d['connect'] for u, v, d in G.edges(data=True)}

    components = list(nx.connected_components(G))
    components_sorted = sorted(components, key=len, reverse=True)

    def _draw_network(graph, layout, colors, labels, edge_label_dict, title):
        fig, ax = plt.subplots(figsize=(12, 8))

        nx.draw_networkx_nodes(
            graph,
            pos=layout,
            node_color=colors,
            cmap=cmap_nodes,
            node_size=node_size,
            alpha=0.8,
            ax=ax
        )
        nx.draw_networkx_edges(
            graph,
            pos=layout,
            edge_color=[a['color'] for _, _, a in graph.edges(data=True)],
            edge_cmap=plt.cm.Set3,
            edge_vmin=0,
            edge_vmax=max(1, len(diffs) - 1),
            width=2,
            ax=ax
        )
        nx.draw_networkx_labels(
            graph,
            pos=layout,
            labels=labels,
            font_size=font_size,
            alpha=0.5,
            ax=ax
        )
        nx.draw_networkx_edge_labels(
            graph,
            pos=layout,
            edge_labels=edge_label_dict,
            font_size=font_size,
            alpha=0.5,
            ax=ax
        )

        finite_colors = pd.Series(colors).dropna().astype(float)
        if len(finite_colors) > 0:
            sm = plt.cm.ScalarMappable(
                cmap=cmap_nodes,
                norm=plt.Normalize(
                    vmin=finite_colors.min(),
                    vmax=finite_colors.max()
                )
            )
            fig.colorbar(sm, ax=ax, label=node_color_label)

        cmap = matplotlib.colormaps['Set3']
        lines = []
        legend_labels = []
        for i, diff in enumerate(diffs):
            rgba = cmap(i / max(1, len(diffs) - 1))
            hex_color = matplotlib.colors.rgb2hex(rgba)
            lines.append(Line2D([0], [0], color=hex_color, linewidth=5))
            legend_labels.append(str(diff))
        ax.legend(lines, legend_labels)

        ax.set_title(title)
        plt.tight_layout()
        plt.show()

    # Plot full network
    n_total = len(mz_values)
    n_in_network = len(G.nodes)
    n_groups = len(components)
    largest_group_size = len(components_sorted[0])
    title_summary = (
        f"n = {n_total}, n_in_network = {n_in_network}, "
        f"n_groups = {n_groups}, largest group = {largest_group_size} nodes"
    )

    fig, ax = plt.subplots(figsize=(12, 8))
    nx.draw_networkx_nodes(
        G,
        pos=pos,
        node_color=graph_node_colors,
        cmap=cmap_nodes,
        node_size=node_size,
        alpha=0.8,
        ax=ax
    )
    nx.draw_networkx_edges(
        G,
        pos=pos,
        edge_color=edge_colors,
        edge_cmap=plt.cm.Set3,
        edge_vmin=0,
        edge_vmax=max(1, len(diffs) - 1),
        width=2,
        ax=ax
    )
    nx.draw_networkx_labels(
        G,
        pos=pos,
        labels=graph_node_labels,
        font_size=font_size,
        alpha=0.5,
        ax=ax
    )
    nx.draw_networkx_edge_labels(
        G,
        pos=pos,
        edge_labels=edge_labels,
        font_size=font_size,
        alpha=0.5,
        ax=ax
    )

    finite_colors = pd.Series(graph_node_colors).dropna().astype(float)
    if len(finite_colors) > 0:
        sm = plt.cm.ScalarMappable(
            cmap=cmap_nodes,
            norm=plt.Normalize(vmin=finite_colors.min(), vmax=finite_colors.max())
        )
        fig.colorbar(sm, ax=ax, label=node_color_label)

    cmap = matplotlib.colormaps['Set3']
    lines = []
    labels = []
    for i, diff in enumerate(diffs):
        rgba = cmap(i / max(1, len(diffs) - 1))
        hex_color = matplotlib.colors.rgb2hex(rgba)
        lines.append(Line2D([0], [0], color=hex_color, linewidth=5))
        labels.append(str(diff))
    ax.legend(lines, labels)

    if show_group_number == True:
        for i, component in enumerate(components_sorted):
            x_coords = [pos[node][0] for node in component]
            y_coords = [pos[node][1] for node in component]
            centroid_x = np.mean(x_coords)
            centroid_y = np.mean(y_coords)
            ax.text(
                centroid_x,
                centroid_y,
                str(i + 1),
                fontsize=30,
                color='gray',
                alpha=0.2,
                ha='center',
                va='center',
                weight='bold'
            )

    ax.set_title(title_summary)
    plt.tight_layout()
    plt.show()

    # # Plot selected network
    # if network_n > len(components_sorted):
    #     print(
    #         f"Warning: group_rank {network_n} exceeds number of groups "
    #         f"({len(components_sorted)}). Plotting largest group instead."
    #     )
    #     network_n = 1

    # selected_component = components_sorted[network_n - 1]
    # G_selected = G.subgraph(selected_component)

    # #pos_selected = nx.spring_layout(G_selected, k=k)
    # pos_selected = gnps_layout(G_selected)
    # node_colors_selected = [
    #     node_color_values[mz_to_idx[node]]
    #     for node in G_selected.nodes()
    # ]
    # node_labels_selected = {
    #     node: graph_node_labels[node]
    #     for node in G_selected.nodes()
    # }
    # edge_labels_selected = {
    #     (u, v): d['connect']
    #     for u, v, d in G_selected.edges(data=True)
    # }

    # _draw_network(
    #     G_selected,
    #     pos_selected,
    #     node_colors_selected,
    #     node_labels_selected,
    #     edge_labels_selected,
    #     f"Network {network_n} has {len(G_selected.nodes())} nodes"
    # )

    # Summary of mass difference detections
    mass_diffs = [data['connect'] for _, _, data in G.edges(data=True)]
    diff_counts = Counter(mass_diffs)

    diffs_detected = [str(d) for d in diff_counts.keys()]
    counts = [diff_counts[d] for d in diff_counts.keys()]

    _, ax_bar = plt.subplots(figsize=(6, 4))
    ax_bar.bar(diffs_detected, counts, color='teal', edgecolor='black', alpha=0.8)
    ax_bar.set_ylabel('Frequency')
    ax_bar.set_xlabel('Mass Difference')
    ax_bar.set_title('Mass Difference Counts in Network')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    return df, diff_matrix

#if __name__ == "__main__":
#    # Example usage
#    mz_array = pd.read_excel('data/mzs_salt_peak_vec_1st.xlsx')['mz'].values
#    mass_difference_network(mz_array, diffs=['H2O', 'NO3', 'HNO3', 'SO4', 'H2SO4', 'HCl', 'Cl'], mz_tol=0.004)


""" 
# Saved code snippets (should be removed if not needed) 

_, idx = np.unique(np.array(edge_colors), return_index=True)
cmap = plt.cm.get_cmap('Set3', len(idx))
c_HEX = []
for i in range(cmap.N):
    rgba = cmap(i)
    c_HEX.append(matplotlib.colors.rgb2hex(rgba))
c_HEX = np.array(c_HEX)[np.array(edge_colors)[np.sort(idx)]]
diffs_new = np.array(diffs)[np.array(edge_colors)[np.sort(idx)]]
lines = [Line2D([0], [0], color=c, linewidth=5) for c in c_HEX]
plt.legend(lines, diffs_new)
plt.show()

G_idxs = nx.from_pandas_edgelist(Df_network, source='idx_1', target='idx_2')
groups = list(nx.connected_components(G_idxs))
groups_list = [list(group) for group in groups]

g = np.array([])
for key in nx.get_edge_attributes(G, 'color').keys():
    g = np.append(g, nx.get_edge_attributes(G, 'color').get(key))

c = []
for edge in G.edges(data=True):
    c.append(edge[2].get('color'))

groups = list(nx.connected_components(G))
for group in groups:
    print(group)

from pyvis.network import Network
import networkx as nx

nx.draw(G, with_labels = True)
nt = Network('500px', '500px')
nt.from_nx(G)
nt.show('nx.html')
"""