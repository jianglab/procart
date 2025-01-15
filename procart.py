""" 
MIT License

Copyright (c) 2021-2024 Wen Jiang

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
required_packages = "streamlit atomium numpy bokeh shapely".split()
import_with_auto_install(required_packages)

import streamlit as st
import numpy as np
np.bool8 = bool  # fix for bokeh 2.4.3
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
        st.elements.lib.policies._shown_default_value_warning = True
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
            fileobj = st.file_uploader(label, type=['pdb', 'cif'], help=None, key="file_upload")
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

        rot_z_auto = round(auto_rotation_angle(model), 1)
        rot_z = st.number_input('Rotation around Z-axis (°)', value=rot_z_auto, min_value=-180.0, max_value=180., step=1.0, key="rot_z")
        rot_x = st.number_input('Rotation around X-axis (°)', value=0.0, min_value=-180.0, max_value=180., step=1.0, key="rot_x")

        show_residue_shape = st.radio('Show residues:', options=["Side chain", "Circle", "Circle (const)", "Blank"], horizontal=True, key="show_residue_shape")
        color_scheme_container = st.container()

        plot_z_dist = st.checkbox('Plot Z-postions of the residues', value=False, key="plot_z_dist")

        with st.expander(label=f"Additional settings", expanded=False):
            save_svg = st.checkbox('Save plots in svg format', value=True, help="save plots in SVG (vector format) if checked or png (pixel format) if unchecked", key="save_svg")
            show_axes = st.checkbox('Show the axes', value=True, key="show_axes")
            show_backbone = st.checkbox('Show backbone of the chains', value=True, key="show_backbone")
            show_ca = st.checkbox('Show Cα of the residues', value=True, key="show_ca")
            show_aa_indices = st.checkbox('Show amino acid indices', value=True, key="show_aa_indices")
            show_gap = st.checkbox('Show gaps in the model', value=True, key="show_gap")
            show_ssbond = st.checkbox('Show disulfide bonds between cysteines', value=True, key="show_ssbond")
            if show_ssbond:
                use_backbone_setting_ssbond = st.checkbox('Use backbone plotting thickness for disulfide bonds', value=True, key="use_backbone_setting_ssbond")
                if not use_backbone_setting_ssbond:
                    ssbond_thickness = int(st.number_input('Disulfide bond line thickness (pixel)', value=3, min_value=0, step=1, key="ssbond_thickness"))
                    ssbond_endpoint_size = int(st.number_input('Disulfide bond marker size (pixel)', value=6, min_value=0, step=1, key="ssbond_endpoint_size"))
                ssbond_color = st.text_input('Disulifde bond line color', placeholder="grey", value="grey", help="example: grey", key="ssbond_color")
                ssbond_endpoint_color = st.text_input('Disulfide bond marker color', placeholder="black", value="black", help="example: grey", key="ssbond_endpoint_color")
            hide_backbone_in_side_chain_shapes = st.checkbox('Hide backbone atoms (CA,C,N,O) when plotting side chain shapes', value=False, key="hide_backbone_in_side_chain_shapes")
            warn_bad_ca_dist = st.checkbox('Warn bad Cα-Cα distances', value=True, key="warn_bad_ca_dist")
            vflip = st.checkbox('Vertically flip the XY-plot', value=False, key="vflip")
            center_xy = st.checkbox('Center the structure in XY plane', value=True, key="center_xy")
            circle_to_background_empty = st.empty()

            center_z = False
            center_z_per_chain = False
            center_zplot_at = ""
            center_zplot_at_aa = []
            one_z_plot = False
            label_at_top = False
            if plot_z_dist:
                one_z_plot = st.checkbox('Plot all Z-plots in one figure', value=True, key="one_z_plot")
                label_at_top = st.checkbox('Label amino acid at the top of Z-plots', value=True, key="label_at_top")
                equal_x = st.checkbox('Use equal spacing in X-axis between amino acids', value=True, key="equal_x")
                center_z_container = st.container()
                center_zplot_at_container = st.container()
                center_zplot_at = center_zplot_at_container.text_input('Center the Z-plot at residues', placeholder="A.1 B.2", value="", help="Center the Z-plot at this amino acid of each chain", key="center_zplot_at")
                if center_zplot_at:
                    center_zplot_at_aa_x = center_z_container.checkbox('Center X-position of the Z-plot', value=True, help="", key="center_zplot_at_aa_x")
                    center_zplot_at_aa_z = center_z_container.checkbox('Center Z-position of the Z-plot', value=True, help="", key="center_zplot_at_aa_z")
                else:
                    center_z = center_z_container.checkbox('Center the structure in Z direction', value=True, key="center_z")

            example = "A: 1-10 17"
            example+= "\nA,B: 1-10 17"
            select_aa = st.text_area("Only plot these residues:", value="", height=72, max_chars=None, key="select_aa", help=None, placeholder=example)

            circle_size_scale = 1.0
            circle_line_thickness = 1
            circle_opaque = 0.9
            if show_residue_shape not in ["Blank"]:
                circle_size_scale = st.number_input('Scale circles relative to the residue sizes', value=1.0, min_value=0.1, step=0.1, key="circle_size_scale")
                circle_line_thickness = int(st.number_input('Circle line width (point)', value=1, min_value=0, step=1, key="circle_line_thickness"))
                circle_opaque = st.number_input('Opaqueness of the circles', value=0.9, min_value=0., max_value=1.0, step=0.1, key="circle_opaque")
                circle_to_background = circle_to_background_empty.checkbox('Send circles to background', value=True, key="circle_to_background")

            aa_label_color = ""
            aa_label_size = 0
            ca_color = ""
            ca_size = 0
            aa_indice_step = 10
            aa_indice_text = ""
            if show_ca:
                ca_size = int(st.number_input('Cα marker size (pixel)', value=6, min_value=0, step=1, key="ca_size"))
                ca_color = st.text_input('Cα color(s)', placeholder="red blue", value="black", help="example: red blue", key="ca_color")
                aa_label_size = int(st.number_input('Amino acid label size (pixel)', value=14, min_value=0, step=1, key="aa_label_size"))
                aa_label_color = st.text_input('Amino acid label color(s)', placeholder="red blue", value="black", help="example: red blue", key="aa_label_color")
                if show_aa_indices:
                    aa_indice_step = int(st.number_input('Show indices every n resiudes', value=10, min_value=1, step=1, key="aa_indice_step"))
                    aa_indice_text = st.text_input("Show indices of these residues:", placeholder="A.13 B.27", value="", help=None, key="aa_indice_text")

            backbone_color = ""
            backbone_thickness = 0
            strand_color = ""
            strand_thickness = 0
            arrowhead_length = 0
            if show_backbone:
                backbone_color = st.text_input('Backbone line color(s)', placeholder="red blue", value="grey", help="example: red blue", key="backbone_color")
                backbone_thickness = int(st.number_input('Backbone line thickness (pixel)', value=3, min_value=0, step=1, key="backbone_thickness"))
                strand_color = st.text_input('Strand line color(s)', placeholder="red blue", value="black", help="example: red blue", key="strand_color")
                strand_thickness = int(st.number_input('Strand line thickness (pixel)', value=6, min_value=0, step=1, key="strand_thickness"))
                arrowhead_length = int(st.number_input('Arrowhead length (pixel)', value=24, min_value=0, step=1, key="arrowhead_length"))
            plot_width = int(st.number_input('Plot width (pixel)', value=1000, min_value=100, step=10, key="plot_width"))
            #transparent_background = st.checkbox('Set background transparent', value=True, key="transparent_background")
            transparent_background = True

        share_url = st.checkbox('Show sharable URL', value=False, help="Include relevant parameters in the browser URL to allow you to share the URL and reproduce the plots", key="share_url")
        if share_url:
            show_qr = st.checkbox('Show QR code of the URL', value=False, help="Display the QR code of the sharable URL", key="show_qr")
        else:
            show_qr = False

        if show_residue_shape != 'Blank' or label_at_top:
            color_scheme = color_scheme_container.radio('Choose a coloring scheme:', options=["Charge", "Hydrophobicity", "Cinema", "Lesk", "Clustal", "Custom"], horizontal=True, key="color_scheme")
            if color_scheme == "Custom":
                example = "white 1-10=red 17=blue L,W=yellow P=cyan"
                example+= "\nA: white 1-10=red 17=blue L,W=yellow P=cyan"
                example+= "\nA,B: white 1-10=red 17=blue L,W=yellow P=cyan"
                custom_color_scheme_txt = color_scheme_container.text_area("Specify your color scheme:", value="", height=128, max_chars=None, key="custom_color_scheme", help=None, placeholder=example)
        else:
            color_scheme = "Cinema"

    chains = []
    for cid in chain_ids:
        chain = model.chain(cid)
        if chain is None:
            st.warning(f"Chain {cid} cannot be found")
        else:
            chains.append((cid, chain.copy()))

    if select_aa:
        chains = select_amino_acids(chains, select_aa)
    
    center_zplot_at_aa = []
    if len(center_zplot_at):
        import re
        center_zplot_at_aa = re.split(r'[;,\s]+', center_zplot_at.upper())
        if len(center_zplot_at_aa) != len(chains):
            center_zplot_at_container.error(f"ERROR: {len(chains)} residues (one per chain) should be specified. You have provided {len(center_zplot_at_aa)} in '{center_zplot_at}'")
            return

    model = atomium.structures.Model(*[chain[1] for chain in chains])
    if rot_z:
        model.rotate(angle=np.deg2rad(rot_z), axis='z')
    if rot_x:
        model.rotate(angle=np.deg2rad(rot_x), axis='x')

    import re
    aa_label_colors = re.split(r'[;,\s]+', aa_label_color)
    ca_colors = re.split(r'[;,\s]+', ca_color)
    backbone_colors = re.split(r'[;,\s]+', backbone_color)
    strand_colors = re.split(r'[;,\s]+', strand_color)

    if vflip: vflip_model(model)

    if center_xy or center_z:
        com_model = model.center_of_mass
        dx = dy = dz = 0.0
        if center_xy:
            dx, dy = -com_model[:2]
        if center_z:
            dz = -com_model[2]
        model.translate(dx=dx, dy=dy, dz=dz)
        if center_z_per_chain:
            for cid, chain in chains:
                com_model = chain.center_of_mass
                chain.translate(dz=-com_model[2])

    from bokeh.plotting import figure
    from bokeh.models import ColumnDataSource, Span, Arrow, VeeHead, HoverTool, Range1d

    fig = figure(x_axis_label="X position (Å)", y_axis_label="Y position (Å)", match_aspect=True)
    if save_svg: fig.output_backend = "svg"
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
    for ci, (cid, chain) in enumerate(chains):
        residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
        res_ids = [f"{res.id}{res.code}" for res in residues]
        seq = [res.code for res in residues]
        ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
        com = np.array([res.center_of_mass for res in residues])
        rog = circle_size_scale*np.array([res.radius_of_gyration for res in residues])
        if show_residue_shape == "Circle (const)":
            rog = np.ones_like(rog) * np.mean(rog)
        strand = [res.strand for res in residues]
        if show_residue_shape  != "Blank":
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

        aa_label_color_i = aa_label_colors[ci%len(aa_label_colors)]
        ca_color_i = ca_colors[ci%len(ca_colors)]
        backbone_color_i = backbone_colors[ci%len(backbone_colors)]
        strand_color_i = strand_colors[ci%len(strand_colors)]

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

        if show_residue_shape in ['Circle', 'Circle (const)']:
            source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos[:,0], 'ca_y':ca_pos[:,1], 'com_x':com[:,0], 'com_y':com[:,1], 'rog':rog, 'color':color, 'strand':strand, 'res_id':res_ids})
            circle = fig.circle(source=source, x='com_x', y='com_y', radius='rog', radius_units='data', line_width=max(1, int(circle_line_thickness)), line_color="black", fill_color='color', fill_alpha=circle_opaque, level='underlay' if circle_to_background else 'overlay')
            hover = HoverTool(renderers=[circle], tooltips=[('COM X', '@com_x{0.00}Å'), ('COM Y', '@com_y{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)
        elif show_residue_shape == 'Side chain':
            def chain_to_bokeh_multi_polygon(chain, scale_factor=1.0, hide_backbone_in_side_chain_shapes=False):
                # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html
                radii = dict(C=1.88, N=1.64, O=1.46, S=1.77, H=1.0, P=1.87, F=1.56, Cl=1.735, Br=1.978, I=2.094)
                if scale_factor!=1.0:
                    for k in radii: radii[k] *= scale_factor
                from shapely.geometry import Point, Polygon
                from shapely.ops import unary_union
                x = []
                y = []
                for res in chain.residues():
                    if not res.atoms(name__regex='CA'): continue
                    if not hide_backbone_in_side_chain_shapes:
                        p = unary_union([Point(a.location[0], a.location[1]).buffer(radii[a.name[0]]) for a in res.atoms()])
                    else:
                        p = unary_union([Point(a.location[0], a.location[1]).buffer(radii[a.name[0]]) for a in set(res.atoms(name__regex='(?!^(CA|O|C|N)$)'))])
                    if not p: continue
                    tmp_x, tmp_y = p.exterior.xy
                    x += [[[tmp_x.tolist()]]]
                    y += [[[tmp_y.tolist()]]]
                return x, y
            mp_x, mp_y = chain_to_bokeh_multi_polygon(chain, scale_factor=circle_size_scale, hide_backbone_in_side_chain_shapes=hide_backbone_in_side_chain_shapes)
            source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos[:,0], 'ca_y':ca_pos[:,1], 'x':mp_x, 'y':mp_y, 'com_x':com[:,0], 'com_y':com[:,1], 'color':color, 'strand':strand, 'res_id':res_ids})
            side_chain = fig.multi_polygons('x', 'y', source=source,  line_width=max(1, int(circle_line_thickness)), line_color="black", fill_color='color', fill_alpha=circle_opaque, level='underlay' if circle_to_background else 'overlay')
            hover = HoverTool(renderers=[side_chain], tooltips=[('COM X', '@com_x{0.00}Å'), ('COM Y', '@com_y{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)

        if backbone_thickness>0:
            source = ColumnDataSource({'x0':nonstrand_x0, 'y0':nonstrand_y0, 'x1':nonstrand_x1, 'y1':nonstrand_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)
            if strand_thickness<=0:
                source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)
                source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)

        if strand_thickness>0:
            source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_thickness, line_color=strand_color_i)

            if arrowhead_length<=0:
                source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_thickness, line_color=strand_color_i)

        if arrowhead_length>0:
            source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
            arrow = Arrow(source=source, x_start='x0', y_start='y0', x_end='x1', y_end='y1', line_width=strand_thickness, line_color=strand_color_i, end=VeeHead(size=arrowhead_length, line_color=strand_color_i, line_width=strand_thickness))
            fig.add_layout(arrow)

        if ca_size>0:
            source = ColumnDataSource({'ca_x':ca_pos[:,0], 'ca_y':ca_pos[:,1], 'res_id':res_ids})
            scatter = fig.scatter(source=source, x='ca_x', y='ca_y', color=ca_color_i, size=ca_size)
            hover = HoverTool(renderers=[scatter], tooltips=[('Cα X', '@ca_x{0.00}Å'), ('Cα Y', '@ca_y{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)

        if show_gap:
            line_color = 'grey'
            source = ColumnDataSource({'x0':gap_x0, 'y0':gap_y0, 'x1':gap_x1, 'y1':gap_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_dash="dotted", line_width=backbone_thickness, line_color=line_color)

        if show_residue_shape != 'Blank' and aa_label_size>0:
            source = ColumnDataSource({'seq':seq, 'com_x':com[:,0], 'com_y':com[:,1]})
            fig.text(source=source, x='com_x', y='com_y', text='seq', text_font_size=f'{aa_label_size:d}pt', text_color="black", text_baseline="middle", text_align="center", level='overlay')

        if show_aa_indices:
            aa_mask = [ ri for ri, res in enumerate(residues) if int(res.id.split('.')[-1])%aa_indice_step==0 ]
            if 0 not in aa_mask: aa_mask = [0] + aa_mask
            if len(residues)-1 not in aa_mask: aa_mask += [len(residues)-1]
            if aa_indice_text:
                import re
                aa_indices_extra = re.split(r'[;,\s]+', aa_indice_text.strip().upper())
                aa_mask += [ri for ri, res in enumerate(residues) if res.id in aa_indices_extra]
            if show_residue_shape != "Blank":
                aa_indices = [residues[i].id.split('.')[-1] for i in aa_mask]
                pos = com
                offset = 0.5*aa_label_size
            else: 
                aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                pos = ca_pos
                offset = 0.5*aa_label_size
            source = ColumnDataSource({'pos_x':pos[aa_mask,0], 'pos_y':pos[aa_mask,1], 'aa_indices':aa_indices})
            fig.text(source=source, x='pos_x', y='pos_y', text='aa_indices', x_offset=offset, text_font_size=f'{aa_label_size:d}pt', text_color=aa_label_color_i, text_baseline="middle", text_align="left", level='overlay')
    
    if show_ssbond:
        if input_mode==0:
            fileobj.seek(0)
            lines = fileobj.getvalue().decode('utf-8').split('\n')
            ssbond_x0 = []
            ssbond_y0 = []
            ssbond_x1 = []
            ssbond_y1 = []
            for line in lines:
                if line.strip()[0:6] == "SSBOND":
                   #print(line)
                   splitted = line.split()
                   if (splitted[3] in chain_ids) and (splitted[6] in chain_ids):
                       atom_0 = model.chain(splitted[3]).residue(splitted[3]+'.'+splitted[4]).atom(name='SG')
                       atom_1 = model.chain(splitted[6]).residue(splitted[6]+'.'+splitted[7]).atom(name='SG')
                       ssbond_x0.append(atom_0.location[0])
                       ssbond_y0.append(atom_0.location[1])
                       ssbond_x1.append(atom_1.location[0])
                       ssbond_y1.append(atom_1.location[1])
            
            if use_backbone_setting_ssbond:
                if show_backbone:
                    final_ssbond_thickness = backbone_thickness
                    final_ssbond_endpoint_size = ca_size                
                else:
                    final_ssbond_thickness = 3
                    final_ssbond_endpoint_size = 6
            else:
                final_ssbond_thickness = ssbond_thickness
                final_ssbond_endpoint_size = ssbond_endpoint_size  

            source = ColumnDataSource({'x0':ssbond_x0, 'y0':ssbond_y0, 'x1':ssbond_x1, 'y1':ssbond_y1})           
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=final_ssbond_thickness, line_color=ssbond_color,line_dash='dashed')   
            
            if ca_size>0:
                scatter = fig.scatter(source=source, x='x0', y='y0', color=ssbond_endpoint_color, size=final_ssbond_endpoint_size)
                scatter = fig.scatter(source=source, x='x1', y='y1', color=ssbond_endpoint_color, size=final_ssbond_endpoint_size)

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
        dxs = []
        for ci, (cid, chain) in enumerate(chains):
            residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
            res_ids = [f"{res.id}{res.code}" for res in residues]
            seq = [res.code for res in residues]
            ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
            com = np.array([res.center_of_mass for res in residues])
            ca_pos_xz, com_xz = unwrap(ca_pos, com, 0 if one_z_plot else center_z) # unwrap the chain to be along x-axis, z-values are preserved
            dxs.append((ca_pos_xz[-1, 0] - ca_pos_xz[0, 0])/(ca_pos_xz.shape[0]-1))
            dxs.append((com_xz[-1, 0]    - com_xz[0, 0])   /(com_xz.shape[0]-1))
            if center_zplot_at_aa:
                center_x_i = None
                for tmp_res_id in range(len(residues)):
                    if residues[tmp_res_id].id == center_zplot_at_aa[ci]:
                        center_x_i = tmp_res_id
                        break
                if center_x_i is not None:
                    if center_zplot_at_aa_x:
                        com_xz[:,0] -= ca_pos_xz[center_x_i][0]
                        ca_pos_xz[:,0] -= ca_pos_xz[center_x_i][0]
                    if center_zplot_at_aa_z:
                        com_xz[:,1] -= ca_pos_xz[center_x_i][1]
                        ca_pos_xz[:,1] -= ca_pos_xz[center_x_i][1]
            ymin = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).min())-1
            ymax = int(np.vstack((ca_pos_xz[:,1], com_xz[:,1])).max())+1
            ymins.append(ymin)
            ymaxs.append(ymax)
        ymin_all = min(ymins)
        ymax_all = max(ymaxs)
        dx_all = np.mean(dxs)

        figs = []
        for ci, (cid, chain) in enumerate(chains):
            residues = [res for res in chain.residues() if res.atoms(name__regex='CA')]
            res_ids = [f"{res.id}{res.code}" for res in residues]
            seq = [res.code for res in residues]
            ca_pos = np.array([list(res.atoms(name__regex='CA'))[0].location for res in residues])
            com = np.array([res.center_of_mass for res in residues])
            ca_pos_xz, com_xz = unwrap(ca_pos, com, 0 if one_z_plot else center_z) # unwrap the chain to be along x-axis, z-values are preserved
            if equal_x:
                ca_pos_xz[:, 0] = ca_pos_xz[0, 0] + dx_all * np.arange(ca_pos_xz.shape[0])
                com_xz[:, 0] = com_xz[0, 0] + dx_all * np.arange(com_xz.shape[0])
            if center_zplot_at_aa:
                center_x_i = None
                for tmp_res_id in range(len(residues)):
                    if residues[tmp_res_id].id == center_zplot_at_aa[ci]:
                        center_x_i = tmp_res_id
                        break
                if center_x_i is not None:
                    if center_zplot_at_aa_x:
                        com_xz[:,0] -= ca_pos_xz[center_x_i][0]
                        ca_pos_xz[:,0] -= ca_pos_xz[center_x_i][0]
                    if center_zplot_at_aa_z:
                        com_xz[:,1] -= ca_pos_xz[center_x_i][1]
                        ca_pos_xz[:,1] -= ca_pos_xz[center_x_i][1]
            rog = circle_size_scale*np.array([res.radius_of_gyration for res in residues])
            strand = [res.strand for res in residues]
            if color_scheme == "Custom":
                res_num = [int(res.id.split(".")[-1]) for res in residues]
                color = np.array(custom_color_mapping(custom_color_scheme_txt, cid, seq, res_num))
            else:
                color = np.array(color_mapping(seq, color_scheme))
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

            aa_label_color_i = aa_label_colors[ci%len(aa_label_colors)]
            ca_color_i = ca_colors[ci%len(ca_colors)]
            backbone_color_i = backbone_colors[ci%len(backbone_colors)]
            strand_color_i = strand_colors[ci%len(strand_colors)]

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
                fig = figure(x_axis_label="Cumulative length of chain projection in XY-plane (Å)", y_axis_label="Cα Z position (Å)", y_range=(ymin, ymax), match_aspect=True)
                figs.append(fig)
                if save_svg: fig.output_backend = "svg"
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
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_dash="dotted", line_width=backbone_thickness, line_color=line_color, level='underlay')

            source = ColumnDataSource({'x0':nonstrand_x0, 'y0':nonstrand_y0, 'x1':nonstrand_x1, 'y1':nonstrand_y1})
            fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)

            if backbone_thickness>0:
                source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
                fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_thickness, line_color=strand_color_i)
                if strand_thickness<=0:
                    source = ColumnDataSource({'x0':strand_body_x0, 'y0':strand_body_y0, 'x1':strand_body_x1, 'y1':strand_body_y1})
                    fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)
                    source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
                    fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=backbone_thickness, line_color=backbone_color_i)

            if strand_thickness>0:
                source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
                arrow = Arrow(source=source, x_start='x0', y_start='y0', x_end='x1', y_end='y1', line_width=strand_thickness, line_color=strand_color_i, end=VeeHead(size=arrowhead_length, line_color=strand_color_i, line_width=strand_thickness))
                fig.add_layout(arrow)

                if arrowhead_length<=0:
                    source = ColumnDataSource({'x0':strand_last_x0, 'y0':strand_last_y0, 'x1':strand_last_x1, 'y1':strand_last_y1})
                    fig.segment(source=source, x0='x0', y0='y0', x1='x1', y1='y1', line_width=strand_thickness, line_color=strand_color_i)

            source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos_xz[:,0], 'ca_z':ca_pos_xz[:,1], 'com_x':com_xz[:,0], 'com_z':com_xz[:,1], 'rog':rog, 'color':color, 'strand':strand, 'res_id':res_ids})
            scatter=fig.scatter(source=source, x='ca_x', y='ca_z', color=ca_color_i, size=ca_size)
            hover = HoverTool(renderers=[scatter], tooltips=[('Chain length', '@ca_x{0.0}Å'), ('Cα Z', '@ca_z{0.00}Å'), ('residue', '@res_id')])
            fig.add_tools(hover)

            if show_residue_shape != "Blank" or label_at_top:
                if label_at_top:
                    if one_z_plot:
                        ca_pos_xz_top = ca_pos_xz[:, 1]*0 + ymax_all + (len(chain_ids) - ci) * (ymax-ymin) * 0.05
                    else:
                        ca_pos_xz_top = ca_pos_xz[:, 1]*0 + ymax + (ymax-ymin) * 0.05
                        fig.y_range=Range1d(ymin, ymax + (ymax-ymin) * 0.05 + 0.5)
                    source = ColumnDataSource({'seq':seq, 'ca_x':ca_pos_xz[:,0], 'ca_z_top':ca_pos_xz_top, 'ca_z':ca_pos_xz[:,1], 'color':color, 'res_id':res_ids})
                    text=fig.text(source=source, x='ca_x', y='ca_z_top', text='seq', x_offset=0, text_font_size=f'{aa_label_size:d}pt', text_color="color", text_baseline="middle", text_align="center")
                    txt_label = f"{residues[0].id.split('.')[-1]}"
                    fig.text(x=[ca_pos_xz[0,0]], y=[ca_pos_xz_top[0]], text=[txt_label], x_offset=-0.5*aa_label_size, y_offset=-0.5*aa_label_size, text_font_size=f'{aa_label_size-4:d}pt', text_color="black", text_baseline="middle", text_align="right", level='overlay')
                    if len(chain_ids)>1: 
                        n_digits = 1+int(np.floor(np.log10(int(txt_label))))
                        txt_label = f"{cid}:"
                        fig.text(x=[ca_pos_xz[0,0]], y=[ca_pos_xz_top[0]], text=[txt_label], x_offset=-(aa_label_size-4)*n_digits, text_font_size=f'{aa_label_size:d}pt', text_color="black", text_baseline="middle", text_align="right", level='overlay')

                    fig.text(x=[ca_pos_xz[-1,0]], y=[ca_pos_xz_top[-1]], text=[f"{residues[-1].id.split('.')[-1]}"], x_offset=0.5*aa_label_size, y_offset=-0.5*aa_label_size, text_font_size=f'{aa_label_size-4:d}pt', text_color="black", text_baseline="middle", text_align="left", level='overlay')
                elif show_residue_shape != "Blank"  and aa_label_size>0:
                    text=fig.text(source=source, x='ca_x', y='ca_z', text='seq', x_offset=aa_label_size, text_font_size=f'{aa_label_size:d}pt', text_color="color", text_baseline="middle", text_align="center", level='overlay')
                    hover = HoverTool(renderers=[text], tooltips=[('Chain length', '@ca_x{0.0}Å'), ('Cα Z', '@ca_z{0.00}Å'), ('residue', '@res_id')])
                    fig.add_tools(hover)
            
            if show_aa_indices:
                aa_mask = [ ri for ri, res in enumerate(residues) if int(res.id.split('.')[-1])%10==0 ]
                if 0 not in aa_mask: aa_mask = [0] + aa_mask
                if len(residues)-1 not in aa_mask: aa_mask += [len(residues)-1]
                if aa_indice_text:
                    import re
                    aa_indices_extra = re.split(r'[;,\s]+', aa_indice_text.strip().upper())
                    aa_mask += [ri for ri, res in enumerate(residues) if res.id in aa_indices_extra]
                if show_residue_shape != "Blank" or label_at_top:
                    if label_at_top:
                        aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                    else:
                        aa_indices = [residues[i].id.split('.')[-1] for i in aa_mask]
                    pos = ca_pos_xz
                    offset = 0.5*aa_label_size if label_at_top else aa_label_size*1.6
                else:
                    aa_indices = [f"{residues[i].code}{residues[i].id.split('.')[-1]}" for i in aa_mask]
                    pos = ca_pos_xz
                    offset = 0.5*aa_label_size
                source = ColumnDataSource({'pos_x':pos[aa_mask,0], 'pos_y':pos[aa_mask,1], 'aa_indices':aa_indices, 'color':color[aa_mask]})
                fig.text(source=source, x='pos_x', y='pos_y', text='aa_indices', x_offset=offset, text_font_size=f'{aa_label_size:d}pt', text_color=aa_label_color_i, text_baseline="middle", text_align="left", level='overlay')

        if len(figs)>1:
            from bokeh.layouts import column
            figs_all = column(children=figs)
            st.bokeh_chart(figs_all)
        elif len(figs)==1:
            if one_z_plot and len(figs):
                figs[0].y_range=Range1d(ymin_all-0.5, ymax_all + (len(chain_ids) * (ymax_all-ymin_all) * 0.05 + 0.5)  * label_at_top)
            st.bokeh_chart(figs[0])

    if share_url:
        set_query_parameters()
        if show_qr:
            qr_image = qr_code()
            st.image(qr_image)
    else:
        st.query_params.clear()

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

def vflip_model(model):
    atoms = model.atoms()
    locations = [[a.location[0], -a.location[1], a.location[2]] for a in atoms]
    for ai, atom in enumerate(atoms):
        atom._location = np.array(locations[ai])

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

def select_amino_acids(chains, select_aa):
    from itertools import product
    keep = set()
    for line in select_aa.split("\n"):
        line = line.strip()
        if len(line)<1: continue
        chain_ids, aa_specs = line.split(":")[:2]   # A,B C: 1-17 23
        chain_ids = [c.upper() for c in chain_ids.replace(",", " ").split()]   # A,B C -> ['A', 'B', 'C']
        indices = []
        for aa_spec in aa_specs.split():  # 1-17 23
            try:
                n1, n2 = aa_spec.split("-")   # 1-17
                indices += list(range(int(n1), int(n2)+1))
            except:
                try:
                    indices += [int(aa_spec)]   # 23
                except:
                    pass
        keep.update(set(map(lambda x: f"{x[0]}.{x[1]}", product(chain_ids, indices))))
    if len(keep)<1: return chains

    ret = []
    for cid, chain in chains:
        residues = [res for res in chain if res.id in keep]
        if residues:
            new_chain = atomium.Chain(*residues)
            ret.append((cid, new_chain))
    return ret

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

# https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
hydrophobicity=dict(I=4.5,V=4.2,L=3.8,F=2.8,C=2.5,M=1.9,A=1.8,G=-0.4,T=-0.7,S=-0.8,W=-0.9,Y=-1.3,P=-1.6,H=-1.6,E=-3.5,Q=-3.5,D=-3.5,N=-3.5,K=-3.9,R=-4.5)
# https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/protein-biology/protein-structural-analysis/amino-acid-reference-chart
pI=dict(A=6,R=10.76,N=5.41,D=2.77,C=5.07,E=3.22,Q=5.65,G=5.97,H=7.59,I=6.02,L=5.98,K=9.74,M=5.74,F=5.48,P=6.3,S=5.68,T=5.6,W=5.89,Y=5.66,V=5.96)
def color_mapping(seq, color_scheme="Cinema"):
    ret = ["white"] * len(seq)
    if color_scheme in ["Hydrophobicity", "Charge"]:
        from bokeh.palettes import Turbo256
        if color_scheme == "Hydrophobicity":
            d = hydrophobicity
            pallet = Turbo256[20:226]   # red (hydrophobic) -> blue (hydrophilic)
        else:
            d = pI
            #pallet = matplotlib colormap: seismic   # red (negative) -> white (neutral) -> blue (positive)
            pallet = Turbo256[::-1][30:236]   # red (negative) -> blue (positive)
        vals = np.array(list(d.values()))
        val_min = vals.min()
        val_max = vals.max()
        for i, aa in enumerate(seq):
            val_int = int(round((d[aa]-val_min)/(val_max-val_min)*(len(pallet)-1)))
            ret[i] = pallet[val_int]
        return ret
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

int_types = dict(aa_indice_step=10, aa_label_size=14, arrowhead_length=24, backbone_thickness=3, ca_size=6, center_xy=1, center_z=1, center_zplot_at_aa_x=0, center_zplot_at_aa_z=0, circle_line_thickness=1, circle_to_background=1, equal_x=1, input_mode=2, label_at_top=1, one_z_plot=1, plot_width=1000, plot_z_dist=0, random_pdb_id=0, save_svg=1, share_url=0, show_aa_indices=1, show_axes=1, show_backbone=1, show_ssbond=1, use_backbone_setting_ssbond=1, ssbond_thickness=3, ssbond_endpoint_size=6, hide_backbone_in_side_chain_shapes=0, show_ca=1, show_gap=1, show_qr=0, strand_thickness=6, transparent_background=1, vflip=0, warn_bad_ca_dist=1)
float_types = dict(circle_opaque=0.9, circle_size_scale=1.0, rot_x=0.0, rot_z=0.0)
other_types = dict(aa_label_color="black", backbone_color="grey", ssbond_color="grey", ssbond_endpoint_color="black", ca_color="black", center_zplot_at="", chain_ids=['A'], color_scheme="Charge", custom_color_scheme="", pdb_id="", select_aa="", show_residue_shape="Side chain", strand_color="black", title="ProCart")
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
    st.query_params.update(d)

def parse_query_parameters():
    query_params = st.query_params
    for attr in query_params:
        if attr == "title":
            st.session_state.title = query_params[attr]
        elif attr == "chain_ids":
            st.session_state.chain_ids = query_params.get_all(attr)
        elif attr in int_types:
            st.session_state[attr] = int(float(query_params[attr]))
        elif attr in float_types:
            st.session_state[attr] = float(query_params[attr])
        else:
            st.session_state[attr] = query_params[attr]
    if "title" not in st.session_state:
        st.session_state.title = "ProCart"

@st.cache_data(show_spinner=False, ttl=24*60*60.) # refresh every day
def get_pdb_ids():
    try:
        #url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
        url ="https://s3.rcsb.org/pub/pdb/derived_data/index/entries.idx"
        ds = np.lib.npyio.DataSource(None)
        with ds.open(url) as fp:
            pdb_ids = [line[:4] for line in fp.readlines()[2:] if len(line) > 4]
    except Exception as e:
        print(e)
        pdb_ids = None
    return pdb_ids

def get_random_pdb_id():
    import random
    pdb_ids_all = get_pdb_ids()
    while True:
        pdb_id = random.choice(pdb_ids_all)
        if is_valid_pdb_id(pdb_id):
            return pdb_id

@st.cache_data(show_spinner=False)
def is_valid_pdb_id(pdb_id):
    if len(pdb_id)!=4: return False
    url = f"https://files.rcsb.org/view/{pdb_id.lower()}.cif"
    ds = np.lib.npyio.DataSource(None)
    if not ds.exists(url):
        return False
    return True

@st.cache_data(show_spinner=False)
def get_model_from_pdb(pdb_id):
    import atomium
    model = atomium.fetch(pdb_id)
    return model

@st.cache_data(show_spinner=False)
def get_model_from_url(url):
    url_final = get_direct_url(url)    # convert cloud drive indirect url to direct url
    ds = np.lib.npyio.DataSource(None)
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

@st.cache_data(persist=True, show_spinner=False)
def setup_anonymous_usage_tracking():
    try:
        import pathlib, stat
        index_file = pathlib.Path(st.__file__).parent / "static/index.html"
        index_file.chmod(stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IROTH)
        txt = index_file.read_text()
        if txt.find("gtag/js?")==-1:
            txt = txt.replace("<head>", '''<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-5Y2WT11MQL"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag('js', new Date());gtag('config', 'G-5Y2WT11MQL');</script>''')
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
        url = "https://procart.streamlit.app/"
    elif host == "heroku":
        url = "https://protein-structure-procart.herokuapp.com/"
    else:
        url = f"http://{host}:8501/"
    import urllib
    params = st.query_params
    d = {k:params[k] for k in params}
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
