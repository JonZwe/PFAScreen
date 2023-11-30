# Functions with plotting routines (e.g., KMD, MD/C-m/C, m/z vs. RT) for PFAScreen
import os
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import altair as alt

# m/z vs. RT plot for MSMS alignment validation
def mz_RT_MSMS(
        Df_FeatureData, 
        Df_MS2RawData, 
        idx_in_features, 
        idx_in_MS2RawData,
        Results_folder
        ):
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Df_FeatureData['RT']/60, y=Df_FeatureData['m/z'],
                        mode='markers',
                        marker=dict(color="Navy", size = 15)))
    fig.add_trace(go.Scatter(x=Df_MS2RawData['RT']/60, y=Df_MS2RawData['m/z'],
                        mode='markers',
                        marker=dict(color="Green", size = 10)))

    fig.add_trace(go.Scatter(x=Df_MS2RawData['RT'][idx_in_MS2RawData]/60, y=Df_MS2RawData['m/z'][idx_in_MS2RawData],
                        mode='markers',
                        marker=dict(color="Orange", size = 8)))
    fig.update_layout(showlegend=False,
                        template='plotly',
                        title={'text': f'{len(idx_in_MS2RawData)} of {len(Df_MS2RawData)} MS2 spectra | {len(np.unique(idx_in_features))} of {len(Df_FeatureData)} Features'},
                        xaxis_title="RT",
                        yaxis_title="m/z",
                        font = dict(size = 20))
    fig.write_html(os.path.join(Results_folder, 'plots/RT_mz_MSMS.html'))

# m/z vs. RT plot
def mz_RT(
        Df_FeatureData,
        Results_folder
        ):
    
    fig = px.scatter(Df_FeatureData, 
                     x = Df_FeatureData['RT']/60, 
                     y = 'm/z',
                     color =  np.log10(Df_FeatureData['m/z intens']), color_continuous_scale = px.colors.sequential.RdPu,
                     hover_name = Df_FeatureData.index, 
                     hover_data=['m/z', 'RT', 'm/z intens', 'min Homologues', 'formulas'])

    fig.update_traces(marker=dict(size=22,
                      line=dict(width=1, color='DarkSlateGrey')),
                      selector=dict(mode='markers'),
                      opacity=0.5)
    
    #fig.add_scatter(x = Df_FeatureData['RT'][Df_FeatureData['FORMULA'].notnull()]/60, 
    #                y = Df_FeatureData['m/z'][Df_FeatureData['FORMULA'].notnull()], 
    #                mode='markers',
    #                marker_size=25, 
    #                marker_symbol='hexagram',
    #                marker=dict(color='black'))
    
    fig.update_layout(xaxis_title="RT (min)", yaxis_title="m/z", font=dict(size=22), showlegend=False)
    #fig.show()
    fig.write_html(os.path.join(Results_folder, 'plots/RT_mz.html'))


# MD/C-m/C plot
def MDC_mC_plot(
        Df_FeatureData,
        Results_folder
        ):
    
    # constants
    m_CF = -8.40596e-05
    m_CHF = -0.0005237
    intercept_CF = 0.0010087
    intercept_CHF = 0.0229902

    x = np.linspace(np.min(Df_FeatureData['m/C']), np.max(Df_FeatureData['m/C']), 100)
    y_CF = m_CF * x + intercept_CF
    y_CHF = m_CHF * x + intercept_CHF

    fig1 = px.scatter(Df_FeatureData, x='m/C', y='MD/C', color = np.log10(Df_FeatureData['m/z intens']),#Df_FeatureData['Score_scaled']
                      color_continuous_scale = px.colors.sequential.Turbo,
                      hover_name = Df_FeatureData.index, 
                      hover_data=['m/z', 'RT', 'm/z intens', 'min Homologues', 'formulas'])
    
    fig2 = px.line(x=x, y=y_CF, markers=False)
    fig3 = px.line(x=x, y=y_CHF, markers=False)
    fig4 = px.line(x = np.ones(200)*50, y = np.linspace(np.min(Df_FeatureData['MD/C']), np.max(Df_FeatureData['MD/C']), 200), markers=False)
    
    fig1.update_traces(marker=dict(size=15,
                      line=dict(width=1, color='DarkSlateGrey')),
                      selector=dict(mode='markers'),
                      opacity=0.5)
    
    fig_all = go.Figure(data=fig1.data + fig2.data + fig3.data + fig4.data)
    fig_all.update_layout(xaxis_title="m/C", yaxis_title="MD/C", font=dict(size=22), showlegend=False)
    fig_all.write_html(os.path.join(Results_folder, 'plots/MDC_mC.html'))

# m/C histogram
def mC_histogram(
        Df_FeatureData,
        Results_folder
        ):
    fig = px.histogram(Df_FeatureData, x="m/C")
    fig.update_layout(xaxis_title="m/C", 
                      yaxis_title="Counts", 
                      font=dict(size=22), 
                      showlegend=False)
    fig.update_xaxes(range=[0, 100])
    fig.write_html(os.path.join(Results_folder, 'plots/mC_histogram.html'))


# KMD vs. m/z plot linked together with m/z vs. RT plot
def KMD_plot(
        Df_FeatureData,
        Results_folder,
        mC_limit = 0
        ):

    Df_FeatureData = Df_FeatureData[['m/z','RT','Unique Homologues','HS Number', 'm/z intens', 'min Homologues', 'm/C', 'KMD', 'formulas']]

    HS_pos = Df_FeatureData[Df_FeatureData['min Homologues'] == True]
    HS_neg = Df_FeatureData[Df_FeatureData['min Homologues'] == False]

    HS_pos = HS_pos[HS_pos['m/C'] > mC_limit]

    # Select features of same homologous series on mouseover
    selection = alt.selection_point(on='mouseover', fields=['HS Number'])

    # Set color of features in the same homologous series
    color = alt.condition(selection,
                        alt.Color('HS Number:N', legend=None, scale=alt.Scale(scheme="set1")),
                        alt.value('lightgrey'))

    # Create layer 1: Figure contains features in homologous series
    HS_features = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='m/z',
        y='KMD',
        color=color,
        tooltip=["m/z","RT","Unique Homologues","HS Number", "m/z intens", "formulas"]
    ).properties(
        width=500,
        height=500  
    ).interactive(
    ).add_params(
        selection
    )

    # Create layer 2: Figure contains features not in homologous series
    HS_negative = alt.Chart(HS_neg).mark_circle(size=100, opacity = 0.01).encode(
        x='m/z',
        y='KMD',
        color=alt.value('lightgray'),
        tooltip=["m/z","RT"],
    )

    # Create layer 3: Activated when a certain homologous series is clicked
    selection2 = alt.selection_point(fields=['HS Number'])
    opacity = alt.condition(selection2, alt.value(0), alt.value(1))


    color_click = alt.condition(selection2,
                        alt.Color('set1:N', legend=None, scale=alt.Scale(scheme="set1")),
                        alt.value('lightgray'))

    HS_click = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='m/z',
        y='KMD',
        color=color_click,
        opacity=opacity,
        tooltip=["m/z","RT","Unique Homologues","HS Number", "m/z intens", 'formulas']
    ).properties(
        width=500,
        height=500  
    ).interactive(
    ).add_params(
        selection2
    )
        
    # Combine all 3 layers for first Figure HS_comb 
    HS_comb = HS_features + HS_negative + HS_click
        

    # Create second figure (RT vs. mz)
    # Create layer 1: Highlight on mouseover
    RT_chart1 = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='RT',
        y='m/z',
        color=color
    ).properties(
        height=500,
        width=500
    ).add_params(
        selection
    )

    # Create layer 2: Activate when clicked
    RT_chart2 = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='RT',
        y='m/z',
        color=color_click,
        opacity=opacity,
        tooltip=["m/z","RT","Unique Homologues","HS Number", "m/z intens", 'formulas']
    ).properties(
        width=500,
        height=500  
    ).interactive(
    ).add_params(
        selection2
    )
        
    # Combine both RT layers
    RT_chart = RT_chart1 + RT_chart2

    # Horizontally concatenate mz vs. KMD and mz vs RT and adjust font size
    HS_concat = alt.hconcat(HS_comb, RT_chart).configure_axis(
        labelFontSize=20,
        titleFontSize=20
    )
    HS_concat.save(os.path.join(Results_folder, 'plots/HS.html'))


# Plotting annotated MS2 spectra (diagnostic fragments and fragment mass differences)
def MS2_spectra_plotter(
        Df_FeatureData,
        idx,
        diffs,
        Results_folder,
        font_size = 16
        ):

    fig = go.Figure()
    for n, peak in enumerate(Df_FeatureData['mz_peaks'][idx]):
        fig.add_trace(
            go.Scatter(x = (Df_FeatureData['mz_peaks'][idx][n], Df_FeatureData['mz_peaks'][idx][n]), y = (Df_FeatureData['intens_peaks'][idx][n], 0), 
            mode = 'lines+text', 
            line = dict(color = "black"),
            text = [np.round(Df_FeatureData['mz_peaks'][idx][n], 4)],
            textposition="top center")
            )
    
    for n, peak in enumerate(Df_FeatureData['mz_peaks_diagnostic'][idx]):
        fig.add_trace(
            go.Scatter(x = (Df_FeatureData['mz_peaks_diagnostic'][idx][n], Df_FeatureData['mz_peaks_diagnostic'][idx][n]), y = (Df_FeatureData['intens_peaks_diagnostic'][idx][n], 0),
            mode = 'lines+text',
            line=dict(color="blue", width = 3),
            text = [Df_FeatureData['formula_diagnostic'][idx][n]],
            textposition="bottom right",
            textfont=dict(color="blue"))
            )
    
    for diff in diffs:
        for n, peak in enumerate(Df_FeatureData[f'mz_{diff}'][idx]):
            fig.add_trace(
                go.Scatter(x = (Df_FeatureData[f'mz_{diff}'][idx][n], Df_FeatureData[f'mz_{diff}'][idx][n]), y = (Df_FeatureData[f'intens_{diff}'][idx][n], 0), 
                mode = 'lines+text', 
                line=dict(color="indianred", width = 3),
                text = [np.round(Df_FeatureData[f'mz_{diff}'][idx][n], 4)],
                textposition="top center",
                textfont=dict(color="indianred"))
                )
            
    if type(Df_FeatureData['frag_idx'][idx]) == np.ndarray:

        u, c = np.unique(Df_FeatureData['frag_idx'][idx], return_counts = True)
        prefix_c = np.zeros(len(Df_FeatureData['frag_idx'][idx]))
        for x in range(len(u)):
            prefix_c[np.where(u[x] == Df_FeatureData['frag_idx'][idx])[0]] = np.arange(0, c[x])
        
        for n, peak in enumerate(Df_FeatureData['frag_idx'][idx]):
            fig.add_trace(
                go.Scatter(x = (Df_FeatureData['mz_peaks'][idx][peak], Df_FeatureData['mz_peaks'][idx][peak]),
                           y = (Df_FeatureData['intens_peaks'][idx][peak], 0),
                mode = 'text', 
                line=dict(color="indianred",width = 3, dash ='dash'),
                text = [int(prefix_c[n])* '<br>' + '<br>' + Df_FeatureData['new_formulas'][idx][n]],
                textposition="bottom right",
                textfont=dict(color="indianred"))
                )
            
    fig.update_layout(showlegend=False,
                      template='plotly', #plotly_white, simple_white
                      title={'text': f'm/z = {np.round(Df_FeatureData["m/z_MSMS"][idx],4)} | RT = {np.round(Df_FeatureData["RT_MSMS"][idx]/60, 2)} | Intensity = {int(Df_FeatureData["intensity"][idx])}'},
                      xaxis_title="m/z",
                      yaxis_title="Counts",
                      xaxis_range=[np.min(Df_FeatureData['mz_peaks'][idx]) - 10, Df_FeatureData["m/z"][idx] + 10],
                      font = dict(size = font_size)
                      )
    fig.write_html(os.path.join(Results_folder, f'plots/Spec_mz_{str(np.round(Df_FeatureData["m/z"][idx], 4))}_intens_{str(np.round(Df_FeatureData["m/z intens"][idx], 0))}.html'))