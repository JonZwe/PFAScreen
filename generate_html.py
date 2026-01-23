
import io
import os
import contextlib
import base64
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from visualization import (isotope_mirror_plot, 
                            spectrum_plot, 
                            mz_rt_plot_interactive, 
                            mdc_mc_plot_interactive, 
                            mc_histogram_plot_interactive, 
                            kmd_plot_interactive)

@contextlib.contextmanager
def suppress_output():
    with open(os.devnull, 'w') as fnull:
        with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
            yield

def sanitize_cell(val):
    try:
        if pd.isna(val):
            return ""
    except:
        pass

    if isinstance(val, float):
        return round(val, 4)
    if isinstance(val, list):
        return [round(x, 4) if isinstance(x, float) else x for x in val]
    return val

def get_columns_for_display(df):
    return [col for col in df.columns if col not in [
        "mzs_isotopes", "ints_isotopes", "mzs_isotopes_theor", "ints_isotopes_theor",
        "mzs_ms2", "ints_ms2", 'SMILES', 'exact_masses', 'adduct_mass_diff', 'annot', 'top_n_idx_ms2'
    ]]

def fig_to_base64(fig):
    with suppress_output():
        buf = io.BytesIO()
        # Reduced DPI for faster generation
        fig.savefig(buf, format='png', dpi=90)
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode('utf-8')

def smiles_to_image_base64(smiles):
    try:
        with suppress_output():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ""
            img = Draw.MolToImage(mol, size=(180, 180))
            buf = io.BytesIO()
            img.save(buf, format='PNG')
            return base64.b64encode(buf.getvalue()).decode('utf-8')
    except:
        return ""

def generate_sample_intensity_plot(row, sample_names):
    try:
        with suppress_output():
            plt.figure(figsize=(3, 2))
            colors = plt.cm.tab10(np.arange(len(sample_names)) % 10)
            values = row[sample_names].values
            plt.bar(range(len(sample_names)), values, color=colors)
            plt.xticks(range(len(sample_names)), sample_names, rotation=90, fontsize=7)
            plt.ylabel('Intensity', fontsize=8)
            plt.tight_layout()
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            plt.close()
            return base64.b64encode(buf.getvalue()).decode('utf-8')
    except:
        return ""

def generate_sample_intensity_plot_fast(row, sample_names):
    """Optimized version for faster HTML generation"""
    try:
        with suppress_output():
            # Use smaller figure and lower DPI for faster generation
            plt.figure(figsize=(2, 1.5))
            colors = plt.cm.tab10(np.arange(len(sample_names)) % 10)
            values = row[sample_names].values
            plt.bar(range(len(sample_names)), values, color=colors)
            plt.xticks(range(len(sample_names)), sample_names, rotation=90, fontsize=6)
            plt.ylabel('Intensity', fontsize=7)
            plt.tight_layout()
            
            buf = io.BytesIO()
            # Reduced DPI for faster generation
            plt.savefig(buf, format='png', dpi=75, bbox_inches='tight')
            plt.close()
            return base64.b64encode(buf.getvalue()).decode('utf-8')
    except:
        return ""

    
def generate_full_html(df, sample_names, output_file="PFAScreen_interactive_table.html", table_name="PFAScreen Table"):

    print('HTML generation started!')

    html_rows = []
    modal_blocks = []
    formula_plot_cache = {}
    display_columns = get_columns_for_display(df)

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Generating HTML"):

        isotope_imgs = []
        smiles_imgs = []
        ms2_preview = ''
        modal_content_ms2 = ''
        sample_intensity_img = ''

        # Generate sample intensity bar plot
        if sample_names is not None and len(sample_names) > 0:
            sample_intensity_base64 = generate_sample_intensity_plot(row, sample_names)
            if sample_intensity_base64:
                sample_intensity_img = f"<img src=\"data:image/png;base64,{sample_intensity_base64}\" style=\"height:80px;cursor:pointer;\" onclick=\"openModal('modal_sample_{idx}')\">"
                # Add modal content for sample intensity plot
                modal_blocks.append(f"<div id=\"modal_sample_{idx}\" class=\"modal\" onclick=\"closeModal('modal_sample_{idx}')\"><img src=\"data:image/png;base64,{sample_intensity_base64}\" style=\"margin:5px;max-width:90vw;max-height:90vh;\"></div>")

        try:
            mz = row['mz']
            mzs_exp = row['mzs_isotopes']
            ints_exp = row['ints_isotopes']
            mzs_theor_all = row['mzs_isotopes_theor']
            ints_theor_all = row['ints_isotopes_theor']
            formulas = row['formulas']
            scores = row['isotope_scores']
            smiles_list = row['SMILES']
            adduct_mass_diff = row['adduct_mass_diff']

            seen = set()
            for n in range(len(formulas)):
                formula = formulas[n]
                delta = adduct_mass_diff[n] if isinstance(adduct_mass_diff, list) else adduct_mass_diff
                key = (formula, delta)

                if key not in seen:
                    if key not in formula_plot_cache:
                        mzs_theor = [x + delta for x in mzs_theor_all[n]]
                        fig = isotope_mirror_plot(mzs_exp, np.array(ints_exp)/max(ints_exp), 
                                                  mzs_theor, ints_theor_all[n], 
                                                  title=f'{formula} | {np.round(scores[n], 3)}')
                        formula_plot_cache[key] = fig_to_base64(fig)
                    isotope_imgs.append(formula_plot_cache[key])
                    seen.add(key)

                smiles_imgs.append(smiles_to_image_base64(smiles_list[n]))
        except:
            pass


        try:
            mz = row['mz']
            mzs_ms2 = row['mzs_ms2']
            ints_ms2 = row['ints_ms2']
            formulas = row['formulas']
            matches_ms2 = row['annot']
            top_n_idx = row['top_n_idx_ms2']

            ms2_imgs = []

            if isinstance(mzs_ms2, list):

                # If matches are NaN or not a list, wrap in one empty dict
                if not isinstance(matches_ms2, list):
                    matches_ms2 = [{}]

                for i, match_dict in enumerate(matches_ms2):
                    # Even if match_dict is {}, plot will still render unannotated spectrum
                    if not isinstance(match_dict, dict):
                        match_dict = {}

                    if isinstance(formulas, list) and i < len(formulas):
                        formula_label = formulas[i]
                    else:
                        formula_label = None

                    fig = spectrum_plot(
                        np.array(mzs_ms2),
                        np.array(ints_ms2),
                        matches=match_dict,
                        top_n_idx=top_n_idx if isinstance(top_n_idx, list) else list(range(len(mzs_ms2))),
                        title=f"MS2 of m/z = {mz} | {formula_label}"
                    )
                    ms2_img = fig_to_base64(fig)
                    ms2_imgs.append(ms2_img)


            # Only show first one as preview
            if ms2_imgs:
                ms2_preview = f"<img src=\"data:image/png;base64,{ms2_imgs[0]}\" style=\"height:100px;cursor:pointer;\" onclick=\"openModal('modal_ms2_{idx}')\">"
                modal_content_ms2 = ''.join(f'<img src="data:image/png;base64,{img}" style="margin:5px;max-width:90vw;max-height:90vh;"><br>' for img in ms2_imgs)
                modal_blocks.append(f"<div id=\"modal_ms2_{idx}\" class=\"modal\" onclick=\"closeModal('modal_ms2_{idx}')\">{modal_content_ms2}</div>")

        except Exception as e:
            print(f"MS2 plot error at row {idx}: {e}")


        isotope_preview = f"<img src=\"data:image/png;base64,{isotope_imgs[0]}\" style=\"height:100px;cursor:pointer;\" onclick=\"openModal('modal_iso_{idx}')\">" if isotope_imgs else ''
        structure_preview = f"<img src=\"data:image/png;base64,{smiles_imgs[0]}\" style=\"height:100px;cursor:pointer;\" onclick=\"openModal('modal_smi_{idx}')\">" if smiles_imgs else ''

        modal_content_iso = ''.join(f'<img src="data:image/png;base64,{img}" style="margin:5px;max-width:90vw;max-height:90vh;"><br>' for img in isotope_imgs)
        modal_content_smi = ''.join(f'<img src="data:image/png;base64,{img}" style="margin:5px;max-width:90vw;max-height:90vh;"><br>' for img in smiles_imgs)

        modal_blocks.append(f"<div id=\"modal_iso_{idx}\" class=\"modal\" onclick=\"closeModal('modal_iso_{idx}')\">{modal_content_iso}</div>")
        modal_blocks.append(f"<div id=\"modal_smi_{idx}\" class=\"modal\" onclick=\"closeModal('modal_smi_{idx}')\">{modal_content_smi}</div>")

        all_cells = ''.join(f"<td>{sanitize_cell(row[col])}</td>" for col in display_columns)

        html_rows.append(f"""
        <tr>
            <td style=\"vertical-align: middle;\">{idx}</td>
            <td style=\"vertical-align: middle;\">{sample_intensity_img}</td>
            <td style=\"vertical-align: middle;\">{isotope_preview}</td>
            <td style=\"vertical-align: middle;\">{structure_preview}</td>
            <td style=\"vertical-align: middle;\">{ms2_preview}</td>
            {all_cells}
        </tr>
        """)

    # Generate plots for the tabs
    mz_rt_html = mz_rt_plot_interactive(df)
    mdc_mc_html = mdc_mc_plot_interactive(df)
    mc_histo = mc_histogram_plot_interactive(df)
    kmd = kmd_plot_interactive(df)

    table_header = ("<th onclick='sortTable(0)'>Index <span class='sort-arrows'>&#9650;&#9660;</span></th><th>Sample Intensities</th><th>Isotopes</th><th>Structures</th><th>MS2</th>" +
                    ''.join(f"<th onclick='sortTable({i + 5})'>{col} <span class='sort-arrows'>&#9650;&#9660;</span></th>"
                            for i, col in enumerate(display_columns)))

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{table_name}</title>
        <style>
            * {{ box-sizing: border-box; }}
            html, body {{
                margin: 0;
                padding: 0;
                height: 100%;
                font-family: Arial, sans-serif;
                overflow: hidden;
            }}
            body > .tab-content {{
                height: calc(100vh - 110px);
                display: flex;
                flex-direction: column;
            }}
            .table-container {{
                flex: 1;
                overflow: auto;
                min-height: 0;
            }}
            table {{ border-collapse: collapse; min-width: max-content; }}
            th, td {{ border: 1px solid black; padding: 5px; text-align: center; vertical-align: middle; }}
            th {{ background-color: #f2f2f2; position: sticky; top: 0; cursor: pointer; }}
            img {{ max-height: 140px; }}
            .sort-arrows {{ font-size: 12px; margin-left: 5px; }}
            .tab-content {{ padding: 20px; }}
            button {{ font-size: 16px; margin-right: 10px; padding: 8px 12px; }}
            .modal {{
                display:none; position:fixed; z-index:1000; left:0; top:0; width:100%; height:100%;
                overflow:auto; background-color:rgba(0,0,0,0.8); text-align:center; padding-top:60px;
            }}
            .modal img {{ max-width:90vw; max-height:90vh; margin:10px auto; }}
        </style>
        <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
        <script>
            let sortDirection = true;
            function sortTable(colIndex) {{
                const table = document.querySelector('table');
                const rows = Array.from(table.rows).slice(1);
                rows.sort((a, b) => {{
                    let A = a.cells[colIndex].innerText;
                    let B = b.cells[colIndex].innerText;
                    let numA = parseFloat(A); let numB = parseFloat(B);
                    if (!isNaN(numA) && !isNaN(numB)) {{ A = numA; B = numB; }}
                    return sortDirection ? A > B ? 1 : -1 : A < B ? 1 : -1;
                }});
                rows.forEach(row => table.appendChild(row));
                sortDirection = !sortDirection;
            }}
            function openModal(id) {{ document.getElementById(id).style.display = 'block'; }}
            function closeModal(id) {{ document.getElementById(id).style.display = 'none'; }}
            function showTab(tabId) {{
                document.querySelectorAll('.tab-content').forEach(el => el.style.display = 'none');
                document.getElementById(tabId).style.display = 'block';
                
                // Trigger Plotly resize when tab becomes visible
                setTimeout(() => {{
                    if (typeof Plotly !== 'undefined') {{
                        Plotly.Plots.resize();
                    }}
                }}, 100);
            }}
            
            // Ensure plots are properly sized when page loads
            window.addEventListener('load', function() {{
                // Check if Plotly loaded successfully
                if (typeof Plotly === 'undefined') {{
                    if (typeof console !== 'undefined') {{
                        console.error('Plotly failed to load from CDN');
                    }}
                    return;
                }}
                
                setTimeout(() => {{
                    Plotly.Plots.resize();
                }}, 500);
            }});
        </script>
    </head>
    <body>

        <h1 style="font-family: Arial; padding: 10px 18px 4px 18px; margin-bottom: 4px;">PF&#916;Screen Report: {table_name}</h1>

        <!-- Tab buttons -->
        <div style="padding: 10px;">
            <button onclick="showTab('table-tab')">Table</button>
            <button onclick="showTab('mzrt-tab')">m/z vs. RT</button>
            <button onclick="showTab('mdcmc-tab')">MD/C vs. m/C</button>
            <button onclick="showTab('mc_histo-tab')">m/C histogram</button>
            <button onclick="showTab('kmd-tab')">KMD</button>
        </div>

        <!-- Table tab -->
        <div id="table-tab" class="tab-content">
            <div class="table-container">
                <table>
                    <tr>{table_header}</tr>
                    {''.join(html_rows)}
                </table>
            </div>
        </div>

        <div id="mzrt-tab" class="tab-content" style="display:none;">
            <h2 style="font-family:Arial">m/z vs. RT</h2>
            {mz_rt_html}
        </div>

        <div id="mdcmc-tab" class="tab-content" style="display:none;">
            <h2 style="font-family:Arial">MD/C vs. m/C</h2>
            {mdc_mc_html}
        </div>

        <div id="mc_histo-tab" class="tab-content" style="display:none;">
            <h2 style="font-family:Arial">m/C histogram</h2>
            {mc_histo}
        </div>

        <div id="kmd-tab" class="tab-content" style="display:none;">
            <h2 style="font-family:Arial">KMD</h2>
            {kmd}
        </div>

        <!-- Modals for zoomed-in views -->
        {''.join(modal_blocks)}

    </body>
    </html>
    """

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"HTML file saved as: {output_file}")