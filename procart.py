""" 
MIT License

Copyright (c) 2021-2022 Wen Jiang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

def import_with_auto_install(packages, scope=locals()):
    if isinstance(packages, str): packages=[packages]
    for package in packages:
        if package.find(":")!=-1:
            package_import_name, package_pip_name = package.split(":")
        else:
            package_import_name, package_pip_name = package, package
        try:
            scope[package_import_name] = __import__(package_import_name)
        except ImportError:
            import subprocess
            subprocess.call(f'pip install {package_pip_name}', shell=True)
            scope[package_import_name] =  __import__(package_import_name)
required_packages = "streamlit atomium numpy bokeh".split()
import_with_auto_install(required_packages)

import streamlit as st
import numpy as np
import atomium

#from memory_profiler import profile
#@profile(precision=4)
def main():
    title = "ProCart"
    st.set_page_config(page_title=title, layout="wide")

    hosted, host = is_hosted(return_host=True)
    if hosted and host in ['heroku']:
        st.error(f"This app hosted on Heroku will be unavailable starting November 28, 2022 [when Heroku discontinues free hosting service](https://blog.heroku.com/next-chapter). Please switch to [the same app hosted elsewhere](https://jianglab-procart-streamlit-app-vxbkzh.streamlitapp.com)")

    session_state = st.session_state
    if "input_mode" not in session_state:  # only run once at the start of the session
        st.elements.utils._shown_default_value_warning = True
        parse_query_parameters()
    st.title(st.session_state.title)

    with st.sidebar:
        input_modes = {0:"upload", 1:"url", 2:"PDB ID"}
        help = None
        input_mode = st.radio(label="How to obtain the input PDB file:", options=list(input_modes.keys()), format_func=lambda i:input_modes[i], index=2, horizontal=True, help=help, key="input_mode")
        pdb_ids_all = get_pdb_ids()
        pdb_ids_amyloid = pdb_ids_all
        pdb = None
        if input_mode == 0: # "upload a PDB file":
            label = "Upload a PDB file"
            fileobj = st.file_uploader(label, type=['pdb'], help=None, key="file_upload")
            if fileobj is not None:
                pdb = get_model_from_uploaded_file(fileobj)
        elif input_mode == 1: # "url":
            url_default = "https://files.rcsb.org/download/7MKH.pdb"
            help = "An online url (http:// or ftp://) or a local file path (/path/to/your/structure.pdb)"
            url = st.text_input(label="Input the url of a PDB file:", value=url_default, help=help, key="url")
            with st.spinner(f'Downloading {url.strip()}'):
                pdb = get_model_from_url(url.strip())
        elif input_mode == 2:
            if not pdb_ids_all:
                st.warning("failed to obtained a list of PDB IDs")
                return
            pdb_id_default = "7MKH"
            do_random_pdb_id = st.checkbox("Choose a random PDB ID", value=False, key="random_pdb_id")
            if do_random_pdb_id:
                help = "Randomly select another structure in PDB"
                button_clicked = st.button(label="Change PDB ID", help=help)
                if button_clicked:
                    st.session_state.pdb_id = get_random_pdb_id()
            else:
                help = None
                label = "Input an PDB ID:"
                pdb_id = st.text_input(label=label, value=pdb_id_default, key='pdb_id', help=help)
                pdb_id = pdb_id.upper()
                if pdb_id not in pdb_ids_all:
                    msg = f"{pdb_id} is not a valid PDB entry. Please input a valid id (for example, '{pdb_id_default}') from the {len(pdb_ids_all):,} entries"
                    st.warning(msg)
                    return
                if pdb_id not in pdb_ids_amyloid:
                    msg= f"{pdb_id} is in PDB but not annotated as an amyloid structure" 
                    st.warning(msg)
            if 'pdb_id' in st.session_state: pdb_id = st.session_state.pdb_id
            else: pdb_id = pdb_id_default
            if not is_valid_pdb_id(pdb_id):
                st.warning(f"{pdb_id} is not a valid PDB entry")
                return
            msg = f'PDB: [{pdb_id}](https://www.rcsb.org/structure/{pdb_id})'
            st.markdown(msg)
            with st.spinner(f'Downloading {pdb_id}'):
                pdb = get_model_from_pdb(pdb_id)
            if pdb is None:
                st.warning(f"Failed to download {pdb_id}](https://www.rcsb.org/structure/{pdb_id})")
                return

        if pdb is None: return

        model = pdb.model
        valid_chain_ids = sorted([chain.id for chain in model.chains()])
        if len(valid_chain_ids)<1:
            st.warning(f"No protein chain in the structure")
            return

        if len(valid_chain_ids)>1:
            chain_ids = st.multiselect('Choose one or more chains:', options=valid_chain_ids, default=[valid_chain_ids[0]], key="chain_ids")
        else:
            chain_ids = valid_chain_ids

        if len(chain_ids)<1:
            st.warning("Please select at least one chain")
            return

        chains = []
        for cid in chain_ids:
            chain = model.chain(cid)
            if chain is None:
                st.warning(f"Chain {cid} cannot be found")
            else:
                chains.append((cid, chain.copy()))
        
        model = atomium.structures.Model(*[chain[1] for chain in chains])
        rotz_auto = round(auto_rotation_angle(model), 1)
        rotz = st.number_input('Rotation around Z-axis (°)', value=rotz_auto, min_value=-180.0, max_value=180., step=1.0, key="rotz")
        if rotz:
            model.rotate(angle=np.deg2rad(rotz), axis='z')

        show_residue_circles = st.checkbox('Show residues in circles', value=True, key="show_residue_circles")
        plot_z_dist = st.checkbox('Plot Z-postions of the residues', value=False, key="plot_z_dist")

        if show_residue_circles:
            color_scheme = st.radio('Choose a coloring scheme:', options=["Cinema", "Lesk", "Clustal", "Custom"], horizontal=True, key="color_scheme")
            if color_scheme == "Custom":
                example = "white 1-10=red 17=blue L,W=yellow P=cyan"
                example+= "\nA: white 1-10=red 17=blue L,W=yellow P=cyan"
                example+= "\nA,B: white 1-10=red 17=blue L,W=yellow P=cyan"
                custom_color_scheme_txt = st.text_area("Specify your color scheme:", value="", height=128, max_chars=None, key="custom_color_scheme", help=None, placeholder=example)
        else:
            color_scheme = "Cinema"

        with st.expander(label=f"Additional settings", expanded=False):
            show_aa_indices = st.checkbox('Show amino acid indices', value=True, key="show_aa_indices")
            show_gap = st.checkbox('Show gaps in the model', value=True, key="show_gap")
            vflip = st.checkbox('Vertically flip the XY-plot', value=False, key="vflip")
            center_xy = st.checkbox('Center the structure in XY plane', value=False, key="center_xy")
            if plot_z_dist:
                center_z = st.checkbox('Center the structure in Z direction', value=False, key="center_z")
                one_z_plot = st.checkbox('Plot all Z-plots in one figure', value=True, key="one_z_plot")
                label_at_top = st.checkbox('Place amino acid labels at the top', value=True, key="label_at_top")
            else:
                center_z = False
                one_z_plot = False
                label_at_top = False
            transparent_background = st.checkbox('Set background transparent', value=True, key="transparent_background")
            show_axes = st.checkbox('Show the axes', value=True, key="show_axes")
            warn_bad_ca_dist = st.checkbox('Warn bad Ca-Ca distances', value=True, key="warn_bad_ca_dist")
            plot_width = int(st.number_input('Plot width (pixel)', value=1000, min_value=100, step=10, key="plot_width"))
            if show_residue_circles:
                circle_size_scale = st.number_input('Scale circles relative to the residue sizes', value=1.0, min_value=0.1, step=0.1, key="circle_size_scale")
                circle_line_thickness = int(st.number_input('Circle line width (point)', value=1, min_value=0, step=1, key="circle_line_thickness"))
                circle_opaque = st.number_input('Opaqueness of the circles', value=0.9, min_value=0., max_value=1.0, step=0.1, key="circle_opaque")
                letter_size = int(st.number_input('Size of the letters (point)', value=16, min_value=1, step=1, key="letter_size"))
            else:
                circle_size_scale = 1.0
                circle_line_thickness = 1
                circle_opaque = 0.9
                letter_size = 10
            backbone_line_thickness = int(st.number_input('Backbone line thickness (pixel)', value=2, min_value=0, step=1, key="backbone_line_thickness"))
            strand_line_thickness = int(st.number_input('Strand line thickness (pixel)', value=4, min_value=0, step=1, key="strand_line_thickness"))

        share_url = st.checkbox('Show sharable URL', value=False, help="Include relevant parameters in the browser URL to allow you to share the URL and reproduce the plots", key="share_url")
        if share_url:
            show_qr = st.checkbox('Show QR code of the URL', value=False, help="Display the QR code of the sharable URL", key="show_qr")
        else:
            show_qr = False

    if center_xy or center_z:
        com_model = model.center_of_mass
        dx = dy = dz = 0.0
        if center_xy:
            dx, dy = -com_model[:2]
        if center_z:
            dz = -com_model[2]
        model.translate(dx=dx, dy=dy, dz=dz)

    from bokeh.plotting import figure
    from bokeh.models import ColumnDataSource, Span, Arrow, VeeHead, HoverTool, Range1d

    fig = figure(x_axis_label="X position (Å)", y_axis_label="Y position (Å)", match_aspect=True)
    fig.frame_width=plot_width
    fig.xgrid.visible = False
    fig.ygrid.visible = False
    if not show_axes:
        fig.xaxis.visible = False
        fig.yaxis.visible = False 
    if transparent_background:
        fig.background_fill_color = None
        fig.border_fill_color = None
        fig.outline_line_color = None

    xmins = []
    xmaxs = []
    ymins = []
    ymaxs = []
    bad_ca_dist = []
    for cid, chain in chains:
        residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
        res_ids = [f"{res.id}{res.code}" for res in residues]
        seq = [res.code for res in residues]
        ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
        com = np.array([res.center_of_mass for res in residues])
        if vflip:
            ca_pos[:, 1] *= -1
            com[:, 1] *= -1
        rog = circle_size_scale*np.array([res.radius_of_gyration for res in residues])
        strand = [res.strand for res in residues]
        if show_residue_circles:
            if color_scheme == "Custom":
                res_num = [int(res.id.split(".")[-1]) for res in residues]
                color = np.array(custom_color_mapping(custom_color_scheme_txt, cid, seq, res_num))
            else:
                color = np.array(color_mapping(seq, color_scheme))
        else:
            color = np.array(['black']*len(seq))

        xmin = int(np.vstack((ca_pos[:,0], com[:,0]-rog)).min())-1
        xmax = int(np.vstack((ca_pos[:,0], com[:,0]+rog)).max())+1
        ymin = int(np.vstack((ca_pos[:,1], com[:,1]-rog)).min())-1
        ymax = int(np.vstack((ca_pos[:,1], com[:,1]+rog)).max())+1
        xmins.append(xmin)
        xmaxs.append(xmax)
        ymins.append(ymin)
        ymaxs.append(ymax)

        strand_body_x0 = []
        strand_body_y0 = []
        strand_body_x1 = []
        strand_body_y1 = []
        strand_last_x0 = []
        strand_last_y0 = []
        strand_last_x1 = []
        strand_last_y1 = []
        nonstrand_x0 = []
        nonstrand_y0 = []
        nonstrand_x1 = []
        nonstrand_y1 = []
        gap_x0 = []
        gap_y0 = []
        gap_x1 = []
        gap_y1 = []

        for i in range(len(strand)):
            ca_dist = np.linalg.norm(ca_pos[i]-ca_pos[i-1])
            if 3.5<ca_dist<4.1:
                if strand[i]:
                    if i==len(strand)-1 or (i<len(strand)-1 and not strand[i+1]): # end of a strand
                        strand_last_x0.append( ca_pos[i-1,0] )
                        strand_last_y0.append( ca_pos[i-1,1] )
                        strand_last_x1.append( ca_pos[i,0] )
                        strand_last_y1.append( ca_pos[i,1] )
                    else:
                        strand_body_x0.append( ca_pos[i-1,0] )
                        strand_body_y0.append( ca_pos[i-1,1] )
                        strand_body_x1.append( ca_pos[i,0] )
                        strand_body_y1.append( ca_pos[i,1] )
                else:
                    nonstrand_x0.append( ca_pos[i-1,0] )
                    nonstrand_y0.append( ca_pos[i-1,1] )
                    nonstrand_x1.append( ca_pos[i,0] )
                    nonstrand_y1.append( ca_pos[i,1] )
            else:
                index_current = int(residues[i].id.split(".")[-1])
                index_previous = int(residues[i-1].id.split(".")[-1])
                if index_current-index_previous==1:
                    bad_ca_dist.append( (ca_dist, res_ids[i-1], res_ids[i]) )
                if index_current>index_previous:
                    gap_x0.append( ca_pos[i-1,0] )
                    gap_y0.append( ca_pos[i-1,1] )
                    gap_x1.append( ca_pos[i,0] )
                    gap_y1.append( ca_pos[i,1] )

        source = ColumnDataSource({'ca_x':ca_pos[:,0], 'ca_y':ca_pos[:,1], 'res_id':res_ids})
        scatter = fig.scatter(source=source, x='ca_x', y='ca_y')
        hover = HoverTool(renderers=[scatter], tooltips=[('Ca X', '@ca_x{0.00}Å'), ('Ca Y', '@ca_y{0.00}Å'), ('residue', '@res_id')])
        fig.add_tools(hover)

        if show_gap:
            line_color = 'grey'
            source = ColumnDataSource({'x0':gap_x0, 'y0':gap_y0, 'x1':gap_x1, 'y1':gap_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_dash="dotted", line_width=backbone_line_thickness, line_color=line_color, level='underlay')

        line_color = 'grey'
        source = ColumnDataSource({'x0':nonstrand_x0, 'y0':nonstrand_y0, 'x1':nonstrand_x1, 'y1':nonstrand_y1})
        fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_line_thickness, line_color=line_color, level='underlay')

        line_color = 'black'
        source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
        fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_line_thickness, line_color=line_color, level='underlay')

        source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
        arrow = Arrow(source=source, x_start='x0', y_start='y0', x_end='x1', y_end='y1', line_width=strand_line_thickness, line_color=line_color, end=VeeHead(size=strand_line_thickness*3, line_color="black", fill_color="black", line_width=strand_line_thickness), level='underlay')
        fig.add_layout(arrow)

        if show_residue_circles:
            source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos[:,0], 'ca_y':ca_pos[:,1], 'com_x':com[:,0], 'com_y':com[:,1], 'rog':rog, 'color':color, 'strand':strand, 'res_id':res_ids})
            circle=fig.circle(source=source, x='com_x', y='com_y', radius='rog', radius_units='data', line_width=max(1, int(circle_line_thickness)), line_color="black", fill_color='color', fill_alpha=circle_opaque, level='guide')
            hover = HoverTool(renderers=[circle], tooltips=[('COM X', '@com_x{0.00}Å'), ('COM Y', '@com_y{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)
            fig.text(source=source, x='com_x', y='com_y', text='seq', text_font_size=f'{letter_size:d}pt', text_color="black", text_baseline="middle", text_align="center", level='overlay')

        if show_aa_indices:
            aa_mask = [ ri for ri, res in enumerate(residues) if int(res.id.split('.')[-1])%10==0 ]
            if 0 not in aa_mask: aa_mask = [0] + aa_mask
            if len(residues)-1 not in aa_mask: aa_mask += [len(residues)-1]
            if show_residue_circles:
                aa_indices = [residues[i].id.split('.')[-1] for i in aa_mask]
                pos = com
                offset = 0.5*letter_size
            else: 
                aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                pos = ca_pos
                offset = 0.5*letter_size
            source = ColumnDataSource({'pos_x':pos[aa_mask,0], 'pos_y':pos[aa_mask,1], 'aa_indices':aa_indices})
            fig.text(source=source, x='pos_x', y='pos_y', text='aa_indices', x_offset=offset, text_font_size=f'{letter_size:d}pt', text_color="black", text_baseline="middle", text_align="left", level='overlay')

    if warn_bad_ca_dist and len(bad_ca_dist):
        bad_ca_dist.sort(key=lambda x: abs(x[0]-3.8), reverse=True)
        pair = "pairs" if len(bad_ca_dist)>1 else "pair"
        msg = f"Warning: {len(bad_ca_dist)} {pair} of neighboring residues with Cα-Cα distance significantly different from the expected distance (3.8Å):  \n"
        msg += "  \n".join([f"{p[1]} - {p[2]}: {p[0]:.2f}Å (err = {p[0]-3.8:.2f}Å)" for p in bad_ca_dist])
        st.warning(msg)
    
    fig.x_range=Range1d(min(xmins)-5, max(xmaxs)+5)
    fig.y_range=Range1d(min(ymins)-5, max(ymaxs)+5)
    fig.frame_height = round(fig.frame_width * (fig.y_range.end-fig.y_range.start)/(fig.x_range.end-fig.x_range.start))
    st.bokeh_chart(fig, use_container_width=False)

    if plot_z_dist:
        ymins = []
        ymaxs = []
        for cid, chain in chains:
            residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
            res_ids = [f"{res.id}{res.code}" for res in residues]
            seq = [res.code for res in residues]
            ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
            com = np.array([res.center_of_mass for res in residues])
            ca_pos_xz, com_xz = unwrap(ca_pos, com, 0 if one_z_plot else center_z) # unwrap the chain to be along x-axis, z-values are preserved
            ymin = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).min())-1
            ymax = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).max())+1
            ymins.append(ymin)
            ymaxs.append(ymax)
        ymin_all = min(ymins)
        ymax_all = max(ymaxs)

        figs = []
        for ci, (cid, chain) in enumerate(chains):
            residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
            res_ids = [f"{res.id}{res.code}" for res in residues]
            seq = [res.code for res in residues]
            ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
            com = np.array([res.center_of_mass for res in residues])
            ca_pos_xz, com_xz = unwrap(ca_pos, com, 0 if one_z_plot else center_z) # unwrap the chain to be along x-axis, z-values are preserved
            rog = circle_size_scale*np.array([res.radius_of_gyration for res in residues])
            strand = [res.strand for res in residues]
            if show_residue_circles:
                if color_scheme == "Custom":
                    res_num = [int(res.id.split(".")[-1]) for res in residues]
                    color = np.array(custom_color_mapping(custom_color_scheme_txt, cid, seq, res_num))
                else:
                    color = np.array(color_mapping(seq, color_scheme))
            else:
                color = np.array(['black']*len(seq))
            white_mask = np.where(color=="white")
            color[white_mask] = "grey"

            ymin = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).min())-1
            ymax = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).max())+1

            strand_body_x0 = []
            strand_body_y0 = []
            strand_body_x1 = []
            strand_body_y1 = []
            strand_last_x0 = []
            strand_last_y0 = []
            strand_last_x1 = []
            strand_last_y1 = []
            nonstrand_x0 = []
            nonstrand_y0 = []
            nonstrand_x1 = []
            nonstrand_y1 = []
            gap_x0 = []
            gap_y0 = []
            gap_x1 = []
            gap_y1 = []

            for i in range(len(strand)):
                ca_dist = np.linalg.norm(ca_pos[i]-ca_pos[i-1])
                if 3.5<ca_dist<4.1:
                    if strand[i]:
                        if i==len(strand)-1 or (i<len(strand)-1 and not strand[i+1]): # end of a strand
                            strand_last_x0.append( ca_pos_xz[i-1,0] )
                            strand_last_y0.append( ca_pos_xz[i-1,1] )
                            strand_last_x1.append( ca_pos_xz[i,0] )
                            strand_last_y1.append( ca_pos_xz[i,1] )
                        else:
                            strand_body_x0.append( ca_pos_xz[i-1,0] )
                            strand_body_y0.append( ca_pos_xz[i-1,1] )
                            strand_body_x1.append( ca_pos_xz[i,0] )
                            strand_body_y1.append( ca_pos_xz[i,1] )
                    else:
                        nonstrand_x0.append( ca_pos_xz[i-1,0] )
                        nonstrand_y0.append( ca_pos_xz[i-1,1] )
                        nonstrand_x1.append( ca_pos_xz[i,0] )
                        nonstrand_y1.append( ca_pos_xz[i,1] )
                else:
                    index_current = int(residues[i].id.split(".")[-1])
                    index_previous = int(residues[i-1].id.split(".")[-1])
                    if index_current>index_previous:
                        gap_x0.append( ca_pos_xz[i-1,0] )
                        gap_y0.append( ca_pos_xz[i-1,1] )
                        gap_x1.append( ca_pos_xz[i,0] )
                        gap_y1.append( ca_pos_xz[i,1] )

            if one_z_plot and len(figs):
                fig = figs[-1]
            else:
                fig = figure(x_axis_label="Chain length (Å)", y_axis_label="Ca Z position (Å)", y_range=(ymin, ymax), match_aspect=True)
                figs.append(fig)
                fig.frame_width=plot_width
                fig.xgrid.visible = False
                fig.ygrid.visible = False
                if not show_axes:
                    fig.xaxis.visible = False
                    fig.yaxis.visible = False 
                if transparent_background:
                    fig.background_fill_color = None
                fig.border_fill_color = None
                fig.outline_line_color = None

                hline = Span(location=0, dimension='width', line_dash='dashed', line_color='black', line_width=1)
                fig.add_layout(hline)

            if show_gap:
                line_color = 'grey'
                source = ColumnDataSource({'x0':gap_x0, 'y0':gap_y0, 'x1':gap_x1, 'y1':gap_y1})
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_dash="dotted", line_width=backbone_line_thickness, line_color=line_color, level='underlay')

            line_color = 'grey'
            source = ColumnDataSource({'x0':nonstrand_x0, 'y0':nonstrand_y0, 'x1':nonstrand_x1, 'y1':nonstrand_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_line_thickness, line_color=line_color, level='underlay')

            line_color = 'black'
            source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_line_thickness, line_color=line_color, level='underlay')

            source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
            arrow = Arrow(source=source, x_start='x0', y_start='y0', x_end='x1', y_end='y1', line_width=strand_line_thickness, line_color=line_color, end=VeeHead(size=strand_line_thickness*3, line_color="black", fill_color="black", line_width=strand_line_thickness), level='underlay')
            fig.add_layout(arrow)

            source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos_xz[:,0], 'ca_z':ca_pos_xz[:,1], 'com_x':com_xz[:,0], 'com_z':com_xz[:,1], 'rog':rog, 'color':color, 'strand':strand, 'res_id':res_ids})
            scatter=fig.scatter(source=source, x='ca_x', y='ca_z')
            hover = HoverTool(renderers=[scatter], tooltips=[('Chain length', '@ca_x{0.0}Å'), ('Ca Z', '@ca_z{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)

            if show_residue_circles or label_at_top:
                if label_at_top:
                    if one_z_plot:
                        ca_pos_xz_top = ca_pos_xz[:, 1]*0 + ymax_all + (len(chain_ids) - ci) * (ymax-ymin) * 0.05
                    else:
                        ca_pos_xz_top = ca_pos_xz[:, 1]*0 + ymax + (ymax-ymin) * 0.05
                        fig.y_range=Range1d(ymin, ymax + (ymax-ymin) * 0.05 + 0.5)
                    source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos_xz[:,0], 'ca_z_top':ca_pos_xz_top, 'ca_z':ca_pos_xz[:,1], 'color':color, 'res_id':res_ids})
                    text=fig.text(source=source, x='ca_x', y='ca_z_top', text='seq', x_offset=0, text_font_size=f'{letter_size:d}pt', text_color="color", text_baseline="middle", text_align="center")
                    txt_label = f"{residues[0].id.split('.')[-1]}"
                    fig.text(x=[ca_pos_xz[0,0]], y=[ca_pos_xz_top[0]], text=[txt_label], x_offset=-0.5*letter_size, y_offset=-0.5*letter_size, text_font_size=f'{letter_size-2:d}pt', text_color="black", text_baseline="middle", text_align="right", level='overlay')
                    if len(chain_ids)>1: 
                        n_txt_label = len(txt_label)
                        txt_label = f"{cid}: "
                        fig.text(x=[ca_pos_xz[0,0]], y=[ca_pos_xz_top[0]], text=[txt_label], x_offset=-(n_txt_label-0.5)*(letter_size-2), text_font_size=f'{letter_size:d}pt', text_color="black", text_baseline="middle", text_align="right", level='overlay')
                    fig.text(x=[ca_pos_xz[-1,0]], y=[ca_pos_xz_top[-1]], text=[f"{residues[-1].id.split('.')[-1]}"], x_offset=0.5*letter_size, y_offset=-0.5*letter_size, text_font_size=f'{letter_size-2:d}pt', text_color="black", text_baseline="middle", text_align="left", level='overlay')
                elif show_residue_circles:
                    text=fig.text(source=source, x='ca_x', y='ca_z', text='seq', x_offset=letter_size, text_font_size=f'{letter_size:d}pt', text_color="color", text_baseline="middle", text_align="center", level='overlay')
                    hover = HoverTool(renderers=[text], tooltips=[('Chain length', '@ca_x{0.0}Å'), ('Ca Z', '@ca_z{0.00}Å'), ('residue', '@res_id')])
                    fig.add_tools(hover)
            
            if show_aa_indices:
                aa_mask = [ ri for ri, res in enumerate(residues) if int(res.id.split('.')[-1])%10==0 ]
                if 0 not in aa_mask: aa_mask = [0] + aa_mask
                if len(residues)-1 not in aa_mask: aa_mask += [len(residues)-1]
                if show_residue_circles or label_at_top:
                    if label_at_top:
                        aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                    else:
                        aa_indices = [residues[i].id.split('.')[-1] for i in aa_mask]
                    pos = ca_pos_xz
                    offset = 0.5*letter_size if label_at_top else letter_size*1.6
                else:
                    aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                    pos = ca_pos_xz
                    offset = 0.5*letter_size
                source = ColumnDataSource({'pos_x':pos[aa_mask,0], 'pos_y':pos[aa_mask,1], 'aa_indices':aa_indices, 'color':color[aa_mask]})
                fig.text(source=source, x='pos_x', y='pos_y', text='aa_indices', x_offset=offset, text_font_size=f'{letter_size:d}pt', text_color="color", text_baseline="middle", text_align="left", level='overlay')

        if len(figs)>1:
            from bokeh.layouts import column
            figs_all = column(children=figs)
            st.bokeh_chart(figs_all)
        elif len(figs)==1:
            if one_z_plot and len(figs):
                figs[0].y_range=Range1d(min(ymins), ymax_all + (len(chain_ids) * (ymax_all-ymin_all) * 0.05 + 0.5)  * label_at_top)
            st.bokeh_chart(figs[0])

    if share_url:
        set_query_parameters()
        if show_qr:
            qr_image = qr_code()
            st.image(qr_image)
    else:
        st.experimental_set_query_params()

    st.markdown("*Developed by the [Jiang Lab@Purdue University](https://jiang.bio.purdue.edu/procart). Report problems to Wen Jiang (jiang12 at purdue.edu)*")

    hide_streamlit_style = """
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

def unwrap(ca, com, center_z=False): # unwrap the chain to be along +x axis
    assert(len(ca)==len(com))
    n = len(com)
    ca_xz = np.zeros((n, 2), dtype=float)
    com_xz = np.zeros((n, 2), dtype=float)
    ca_xz[:, 1] = ca[:, 2] # z -> z
    com_xz[:, 1] = com[:, 2] # z -> z
    for i in range(1, n):
        ca_vec = ca[i]-ca[i-1]
        ca_vec_xy_len = np.linalg.norm(ca_vec[:2])
        ca_xz[i, 0] = ca_xz[i-1, 0] +  ca_vec_xy_len # x,y --> x
        ca_vec_xy_norm = ca_vec[:2]/ca_vec_xy_len
        if i==1:
            com_vec = com[i-1]-ca[i-1]
            com_x_proj = com_vec[0]*ca_vec_xy_norm[0] + com_vec[1]*ca_vec_xy_norm[1]
            com_xz[i-1, 0] = ca_xz[i-1, 0] + com_x_proj # x,y --> x
            com_xz[i-1, 1] = ca_xz[i-1, 1] + com_vec[2] # x,y --> x
        com_vec = com[i]-ca[i]
        com_x_proj = com_vec[0]*ca_vec_xy_norm[0] + com_vec[1]*ca_vec_xy_norm[1]
        com_xz[i, 0] = ca_xz[i, 0] + com_x_proj # x,y --> x
    if center_z:
        z_mean = ca_xz[:, 1].mean()
        ca_xz[:, 1] -= z_mean
        com_xz[:, 1] -= z_mean
    return ca_xz, com_xz  

def auto_rotation_angle(model):
    com = model.center_of_mass[:2]
    ca_atoms = model.atoms(name__regex='CA')
    dxy = np.array([atom.location[:2]-com for atom in ca_atoms])
    inertia = np.dot(dxy.transpose(), dxy)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)[::-1]
    e_values = e_values[order]
    e_vectors = e_vectors[:, order].transpose()
    angle = np.rad2deg(np.arctan2(e_vectors[0, 1], e_vectors[0, 0]))
    return -angle

def custom_color_mapping(custom_color_scheme_txt, cid, seq, res_num):
    assert(len(seq) ==  len(res_num))
    ret = ["white"] * len(seq)
    for line in custom_color_scheme_txt.split("\n"):
        line = line.strip()
        if len(line)<1: continue
        try: # when chain ids are provided -> applicable only to those chains
            chains, color_specs = line.split(":")[:2]   # A,B C: white 1-17=red A,L,W=green
            chains = [c.upper() for c in chains.replace(",", " ").split()]   # A,B C -> ['A', 'B', 'C']
        except: # if no chain ids are provided -> applicable to all chains
            chains = None
            color_specs = line
        if chains and cid.upper() not in chains: continue
        for color_spec in color_specs.split():  # white 1-17=red A,L,W=green
            try:
                res_num_letter, color = color_spec.split("=")   # 1-17=red or A,L,W=green
            except:
                res_num_letter = None   # white
                color = color_spec
            if res_num_letter is None:  # white
                ret = [color] * len(seq)
            elif res_num_letter.find("-")!=-1:  # 1-17=red
                try:
                    start, end = [int(num) for num in res_num_letter.split("-")]
                    for i in range(len(seq)):
                        if start <= res_num[i] <= end:
                            ret[i] = color
                except:
                    st.warning(f"ignoring '{color_spec}'\t{res_num_letter}")
            elif res_num_letter.isnumeric():
                res_num_int = int(res_num_letter)
                for i in range(len(seq)):
                    if res_num_int == res_num[i]:
                        ret[i] = color
            else:  # A,L,W=green or A=green
                residues = [l.upper() for l in res_num_letter.split(",")]
                for i in range(len(seq)):
                    if seq[i].upper() in residues:
                        ret[i] = color
    return ret

def color_mapping(seq, color_scheme="Cinema"):
    ret = ["white"] * len(seq)
    for i, aa in enumerate(seq):
        # https://www.bioinformatics.nl/~berndb/aacolour.html
        if color_scheme == "Cinema":
            if aa in 'H K R'.split(): ret[i] = 'cyan'
            elif aa in 'D E'.split(): ret[i] = 'red'
            elif aa in 'S T N Q'.split(): ret[i] = 'green'
            elif aa in 'A V L I M'.split(): ret[i] = 'white'
            elif aa in ['F W Y']: ret[i] = 'magenta'
            elif aa in ['P G']: ret[i] = 'brown'
            elif aa in ['C']: ret[i] = 'yellow'
            else: ret[i] = 'grey'
        elif color_scheme == "Clustal":
            if aa in 'G P S T'.split(): ret[i] = 'orange'
            elif aa in 'H K R'.split(): ret[i] = 'red'
            elif aa in 'F W Y'.split(): ret[i] = 'blue'
            elif aa in 'I L M V'.split(): ret[i] = 'green'
            else: ret[i] = 'white'
        else: # "Lesk":
            if aa in 'G A S T'.split(): ret[i] = 'orange'
            elif aa in 'C V I L P F Y M W'.split(): ret[i] = 'green'
            elif aa in 'N Q H'.split(): ret[i] = 'magenta'
            elif aa in 'D E'.split(): ret[i] = 'red'
            elif aa in 'K R'.split(): ret[i] = 'blue'
            else: ret[i] = 'white'
    return ret

int_types = dict(backbone_line_thickness=2, center_xy=0, center_z=0, circle_line_thickness=1, input_mode=2, label_at_top=1, letter_size=10, one_z_plot=1, plot_width=1000, plot_z_dist=0, random_pdb_id=0, share_url=0, show_aa_indices=1, show_axes=1, show_gap=1, show_qr=0, show_residue_circles=1, strand_line_thickness=4, transparent_background=1, vflip=0, warn_bad_ca_dist=1)
float_types = dict(circle_opaque=0.9, circle_size_scale=1.0, rotz=0.0)
other_types = dict(chain_ids=['A'], color_scheme="Cinema", custom_color_scheme="", title="ProCart")
def set_query_parameters():
    d = {}
    for k in sorted(st.session_state.keys()):
        v = st.session_state[k]
        if k in int_types:
            v=int(v)
            if v==int_types[k]: continue
        elif k in float_types: 
            v=float(v)
            if v==float_types[k]: continue
        elif k in other_types:
            if v==other_types[k]: continue
        d[k] = v
    st.experimental_set_query_params(**d)

def parse_query_parameters():
    query_params = st.experimental_get_query_params()
    for attr in query_params:
        if attr == "title":
            st.session_state.title = query_params[attr][0]
        elif attr == "chain_ids":
            st.session_state.chain_ids = query_params[attr]
        elif attr in int_types:
            st.session_state[attr] = int(float(query_params[attr][0]))
        elif attr in float_types:
            st.session_state[attr] = float(query_params[attr][0])
        else:
            st.session_state[attr] = query_params[attr][0]
    if "title" not in st.session_state:
        st.session_state.title = "ProCart"

@st.cache(persist=True, show_spinner=False, ttl=24*60*60.) # refresh every day
def get_pdb_ids():
    try:
        url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
        ds = np.DataSource(None)
        with ds.open(url) as fp:
            pdb_ids = [line[:4] for line in fp.readlines()[2:] if len(line) > 4]
    except:
        pdb_ids = None
    return pdb_ids

def get_random_pdb_id():
    import random
    pdb_ids_all = get_pdb_ids()
    while True:
        pdb_id = random.choice(pdb_ids_all)
        if is_valid_pdb_id(pdb_id):
            return pdb_id

@st.experimental_singleton(show_spinner=False)
def is_valid_pdb_id(pdb_id):
    if len(pdb_id)!=4: return False
    url = f"https://files.rcsb.org/view/{pdb_id.lower()}.cif"
    ds = np.DataSource(None)
    if not ds.exists(url):
        return False
    return True

@st.experimental_singleton(show_spinner=False)
def get_model_from_pdb(pdb_id):
    import atomium
    model = atomium.fetch(pdb_id)
    return model

@st.experimental_singleton(show_spinner=False, suppress_st_warning=True)
def get_model_from_url(url):
    url_final = get_direct_url(url)    # convert cloud drive indirect url to direct url
    ds = np.DataSource(None)
    if not ds.exists(url_final):
        st.error(f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview")
        st.stop()
    with ds.open(url) as fp:
        model = get_model_from_file(fp.name)
    return model

def get_model_from_file(filename):
    if filename.endswith(".gz"):
        filename_final = filename[:-3]
        import gzip, shutil
        with gzip.open(filename, 'r') as f_in, open(filename_final, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        filename_final = filename
    import atomium
    model = atomium.open(filename_final)
    return model

def get_model_from_uploaded_file(fileobj):
    import os, tempfile
    orignal_filename = fileobj.name
    suffix = os.path.splitext(orignal_filename)[-1]
    with tempfile.NamedTemporaryFile(suffix=suffix) as temp:
        temp.write(fileobj.read())
        return get_model_from_file(temp.name)

def get_direct_url(url):
    import re
    if url.startswith("https://drive.google.com/file/d/"):
        hash = url.split("/")[5]
        return f"https://drive.google.com/uc?export=download&id={hash}"
    elif url.startswith("https://app.box.com/s/"):
        hash = url.split("/")[-1]
        return f"https://app.box.com/shared/static/{hash}"
    elif url.startswith("https://www.dropbox.com"):
        if url.find("dl=1")!=-1: return url
        elif url.find("dl=0")!=-1: return url.replace("dl=0", "dl=1")
        else: return url+"?dl=1"
    elif url.find("sharepoint.com")!=-1 and url.find("guestaccess.aspx")!=-1:
        return url.replace("guestaccess.aspx", "download.aspx")
    elif url.startswith("https://1drv.ms"):
        import base64
        data_bytes64 = base64.b64encode(bytes(url, 'utf-8'))
        data_bytes64_String = data_bytes64.decode('utf-8').replace('/','_').replace('+','-').rstrip("=")
        return f"https://api.onedrive.com/v1.0/shares/u!{data_bytes64_String}/root/content"
    else:
        return url

@st.cache(persist=True, show_spinner=False)
def setup_anonymous_usage_tracking():
    try:
        import pathlib, stat
        index_file = pathlib.Path(st.__file__).parent / "static/index.html"
        index_file.chmod(stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IROTH)
        txt = index_file.read_text()
        if txt.find("gtag/js?")==-1:
            txt = txt.replace("<head>", '''<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-5Y2WT11MQL"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag('js', new Date());gtag('config', 'G-YV3ZFR8VG6');</script>''')
            index_file.write_text(txt)
    except:
        pass

def get_username():
    from getpass import getuser
    return getuser()

def get_hostname():
    import socket
    fqdn = socket.getfqdn()
    return fqdn

def is_hosted(return_host=False):
    hosted = False
    host = ""
    fqdn = get_hostname()
    if fqdn.find("heroku")!=-1:
        hosted = True
        host = "heroku"
    username = get_username()
    if username.find("appuser")!=-1:
        hosted = True
        host = "streamlit"
    if not host:
        host = "localhost"
    if return_host:
        return hosted, host
    else:
        return hosted

def get_url():
    # ad hoc way before streamlit can return the url
    _, host = is_hosted(return_host=True)
    if len(host)<1: return None
    if host == "streamlit":
        url = "https://share.streamlit.io/wjiang/procart/main/"
    elif host == "heroku":
        url = "https://protein-structure-procart.herokuapp.com/"
    else:
        url = f"http://{host}:8501/"
    import urllib
    params = st.experimental_get_query_params()
    d = {k:params[k][0] for k in params}
    url += "?" + urllib.parse.urlencode(d)
    return url

def qr_code(url=None, size = 8):
    import_with_auto_install(["qrcode"])
    import qrcode
    if url is None:
        url = get_url()
    if not url: return None
    img = qrcode.make(url)  # qrcode.image.pil.PilImage
    data = np.array(img.convert("RGBA"))
    return data

def print_memory_usage():
    from inspect import currentframe
    import psutil, os
    cf = currentframe()
    print(f'Line {cf.f_back.f_lineno}: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2} MB')

if __name__ == "__main__":
    setup_anonymous_usage_tracking()
    main()
