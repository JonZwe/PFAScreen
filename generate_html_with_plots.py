import os

def generate_html_with_plots(df, 
                             plot_dir,
                             table_name='EIC-based suspect screening',
                             output_html='interactive_table.html'):
    """
    Generates an HTML file displaying a table with plots.
    Includes isotope plots, MS2 spectrum plots (if available), and structure plots as part of the enlarged view.
    
    Parameters:
    - df (pd.DataFrame): The DataFrame containing data to display in the table.
    - plot_dir (str): Path to the directory containing the plot images.
        Assumes filenames are of the form 'eic_thumb_X.png', 'eic_large_X.png',
        'eic_complete_X.png', 'isotopes_X.png', 'structure_X.png', and optionally 'ms2_spectrum_X.png'.
    - output_html (str): Path to save the generated HTML file.
    """
    # Check if the plot directory exists
    if not os.path.exists(plot_dir):
        raise FileNotFoundError(f"Directory '{plot_dir}' does not exist.")

    # Generate HTML rows
    html_rows = []
    for index, row in df.iterrows():
        # Dynamically generate cells for all columns
        data_cells = "".join(f"<td>{row[col]}</td>" for col in df.columns)

        # Paths for thumbnails, large plots, complete plots, isotopes, structures, and MS2 spectra
        thumbnail_path = os.path.join(plot_dir, f"eic_thumb_{index}.png")
        large_path = os.path.join(plot_dir, f"eic_large_{index}.png")
        complete_path = os.path.join(plot_dir, f"eic_complete_{index}.png")
        isotope_path = os.path.join(plot_dir, f"isotopes_{index}.png")
        ms2_spectrum_path = os.path.join(plot_dir, f"ms2_spectrum_{index}.png")
        structure_path = os.path.join(plot_dir, f"structure_{index}.png")

        # Ensure paths are properly handled for browser compatibility
        thumbnail_path_web = thumbnail_path.replace("\\", "/")
        large_path_web = large_path.replace("\\", "/")
        complete_path_web = complete_path.replace("\\", "/")
        isotope_path_web = isotope_path.replace("\\", "/")
        ms2_spectrum_path_web = ms2_spectrum_path.replace("\\", "/")
        structure_path_web = structure_path.replace("\\", "/")

        # Check if MS2 spectrum exists
        ms2_spectrum_html = f"""
            <img id="ms2-spectrum-plot" src="{ms2_spectrum_path_web}" alt="MS2 Spectrum Plot" 
                 style="margin-top:10px;max-width:75%;max-height:75%;display:block;">
        """ if os.path.exists(ms2_spectrum_path) else ""

        # Add row to HTML
        html_rows.append(f"""
        <tr>
            <td>
                <img src="{thumbnail_path_web}" alt="Plot {index}" 
                     style="cursor:pointer;width:200px;height:auto;"  
                     onclick="showPlots('{large_path_web}', '{complete_path_web}', '{isotope_path_web}', '{ms2_spectrum_path_web}')">
            </td>
            <td>
                <img src="{structure_path_web}" alt="Structure {index}" 
                     style="cursor:pointer;width:100px;height:auto;"  
                     onclick="showStructure('{structure_path_web}')">
            </td>
            {data_cells}
        </tr>
        """)

    # Generate table header
    table_header = "<th>Plot</th><th>Structure</th>" + "".join(f"<th onclick='sortTable({i})'>{col} <span class='sort-arrows'>▲▼</span></th>" for i, col in enumerate(df.columns, start=2))

    # Full HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{table_name}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                table-layout: auto;
            }}
            th, td {{
                border: 1px solid black;
                padding: 5px;
                text-align: center;
                vertical-align: middle;
            }}
            th {{
                background-color: #f2f2f2;
                position: sticky;
                top: 0;
                z-index: 2;
                cursor: pointer;
            }}
            img {{
                display: block;
                margin: auto;
            }}
            .plot-container {{
                display: none;
                position: fixed;
                top: 5%;
                left: 10%;
                width: 80%;
                height: 85%;
                background: white;
                border: 2px solid black;
                box-shadow: 0 0 10px rgba(0, 0, 0, 0.5);
                z-index: 1000;
                overflow: auto;
                text-align: center;
                padding: 20px;
            }}
            .plot-container img {{
                max-width: 75%;
                max-height: 75%;
                margin: auto;
            }}
            .sort-arrows {{
                font-size: 12px;
                margin-left: 5px;
            }}
        </style>
        <script>
            let sortDirection = true; // true for ascending, false for descending
            function sortTable(columnIndex) {{
                const table = document.querySelector("table");
                const rows = Array.from(table.rows).slice(1); // Exclude header
                rows.sort((a, b) => {{
                    const cellA = a.cells[columnIndex].innerText.toLowerCase();
                    const cellB = b.cells[columnIndex].innerText.toLowerCase();
                    if (!isNaN(cellA) && !isNaN(cellB)) {{
                        return sortDirection ? cellA - cellB : cellB - cellA;
                    }}
                    return sortDirection ? cellA.localeCompare(cellB) : cellB.localeCompare(cellA);
                }});
                rows.forEach(row => table.appendChild(row));
                sortDirection = !sortDirection; // Toggle sorting direction
            }}
            function showPlots(largeImagePath, completeImagePath, isotopeImagePath, ms2SpectrumPath) {{
                const container = document.getElementById('plot-container');
                container.style.display = 'block';
                document.getElementById('large-plot').src = largeImagePath;
                document.getElementById('complete-plot').src = completeImagePath;
                document.getElementById('isotope-plot').src = isotopeImagePath;
                const ms2SpectrumPlot = document.getElementById('ms2-spectrum-plot');
                if (ms2SpectrumPath && ms2SpectrumPath !== 'undefined') {{
                    ms2SpectrumPlot.style.display = 'block';
                    ms2SpectrumPlot.src = ms2SpectrumPath;
                }} else {{
                    ms2SpectrumPlot.style.display = 'none';
                }}
            }}
            function showStructure(structureImagePath) {{
                const container = document.getElementById('structure-container');
                container.style.display = 'block';
                document.getElementById('structure-plot').src = structureImagePath;
            }}
            function closePlot() {{
                const container = document.getElementById('plot-container');
                container.style.display = 'none';
            }}
            function closeStructure() {{
                const container = document.getElementById('structure-container');
                container.style.display = 'none';
            }}
        </script>
    </head>
    <body>
        <h1>{table_name}</h1>
        <table>
            <tr>
                {table_header}
            </tr>
            {''.join(html_rows)}
        </table>
        <div id="plot-container" class="plot-container" onclick="closePlot()">
            <img id="large-plot" src="" alt="Large Plot">
            <img id="complete-plot" src="" alt="Complete Plot">
            <img id="isotope-plot" src="" alt="Isotope Plot">
            {ms2_spectrum_html}
        </div>
        <div id="structure-container" class="plot-container" onclick="closeStructure()">
            <img id="structure-plot" src="" alt="Structure Plot">
        </div>
    </body>
    </html>
    """

    # Save the HTML file
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"HTML file saved as {output_html}. Open it in a browser to view.")
