"""Example Dash App."""

import dash
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import ctx, dcc, html
from dash.dependencies import Input, Output, State
from icecream import ic
from monty.serialization import loadfn
from pymatgen.analysis.defects.plotting.thermo import (
    _label_slopes,
    _plot_line,
    get_plot_data,
)

FEDS = [*loadfn("fed_example_plotting_AlN.json"), *loadfn("fed_example_plotting.json")]


def parse_feds(feds: list) -> dict:
    """Process the formation energy data for plotting.

    Given a list of formation energy data, return a dictionary of parsed data.

    Args:
        feds (list): List of formation energy data.

    Returns:
        A dictionary of plot data.
            - key: The unique identifier (just the unique name
                from group_formation_energy_diagrams).
            - value: A dictionary with the following keys:
                - fed: The formation energy diagram.
                - style: The style of the line.
                - color: The color of the line.
                - x_anno: The x-coordinate of the annotation.
                - chempot: The chemical potentials used to generate the transitions.
    """
    red_dict = get_plot_data(feds, chempot=None)
    for key, data_ in red_dict.items():
        red_dict[key]["chemsys"] = data_["fed"].defect_chemsys
    return red_dict


FEDS_DATA = parse_feds(FEDS)
ALL_CHEMSYS = list({data["chemsys"] for data in FEDS_DATA.values()})

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

chemsys_drop = dcc.Dropdown(
    id="chemsys_dropdown",
    options=[{"label": chemsys, "value": chemsys} for chemsys in ALL_CHEMSYS],
    value="Al-N",
)


# draw chempot diagram
@app.callback(Output("chempot", "figure"), Input("chemsys_dropdown", "value"))
def update_chempot_diagram(chemsys: str) -> go.Figure:
    """Cycle through the formation energy data."""
    fed = next(data["fed"] for data in FEDS_DATA.values() if data["chemsys"] == chemsys)
    chempot_diagram = fed.chempot_diagram
    bulk_formula = fed.bulk_formula
    chempot_plot = chempot_diagram.get_plot(
        formulas_to_draw=[bulk_formula],
        draw_formula_meshes=True,
        draw_formula_lines=False,
    )
    chempot_plot.update_layout(showlegend=False)
    if chempot_diagram.dim == 3:
        for d_ in chempot_plot.data:
            if isinstance(d_, go.Mesh3d):
                d_.name = d_.name.removesuffix(" (mesh)")
                continue
            d_.hoverinfo = "skip"
    return chempot_plot


graph_fed = dcc.Graph(id="formation_en")
graph_chempot = dcc.Graph(id="chempot", figure=update_chempot_diagram(ALL_CHEMSYS[0]))

chemsys_card = dbc.Card(
    [
        html.Div(
            [
                dbc.Label("Select chemical system"),
                chemsys_drop,
            ]
        ),
        html.Div([dbc.Label("Y variable"), graph_chempot]),
    ],
    body=True,
)

app.layout = dbc.Container(
    [
        html.H1("TITLE"),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(chemsys_card, md=6),
                dbc.Col(graph_fed, md=4),
            ],
            align="center",
        ),
    ],
    fluid=True,
)


def _plot_fed_data(fed_data: dict, chemsys: str = None) -> go.Figure:
    fig = go.Figure()
    for uid, data in fed_data.items():
        if chemsys is not None and data["chemsys"] != chemsys:
            continue
        _plot_line(
            pts=data["fed"].get_transitions(data["chempot"]),
            fig=fig,
            name=uid,
            color=data["color"],
            style=data["style"],
            meta={"formation_energy_plot": True},
            uid=uid,
        )
    _label_slopes(fig)
    fig.update_layout(
        title="Formation Energy Diagrams",
        xaxis_title="Fermi Level (eV)",
        yaxis_title="Formation Energy (eV)",
        template="plotly_white",
        font_family="Helvetica",
        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=False),
        showlegend=True,
    )
    return fig


def _shift_slope_labels(fig: dict, filter_uids: list[str], shift: float) -> None:
    for data_ in fig["data"]:
        uid_ = data_.get("uid", None)
        if uid_ not in filter_uids:
            continue
        data_["y"] = [y_ + shift for y_ in data_["y"]]


# def _shift_annotation(fig, name, color, shift):
#     for anno in fig["layout"]["annotations"]:
#         if anno["text"] == name and anno["font"]["color"] == color:
#             anno["y"] += shift


def _update_fed_from_2d_chempot(
    fed_data: dict, chemsys: str, d_click_data: dict, formation_en_fig: dict
) -> None:
    """Update all materials with the same chemsys."""
    update_uids = {uid for uid, data in fed_data.items() if data["chemsys"] == chemsys}
    x, y = d_click_data["points"][0]["x"], d_click_data["points"][0]["y"]
    for data_ in formation_en_fig["data"]:
        uid_ = data_.get("uid", None)
        if uid_ not in update_uids:
            continue
        fed_ = fed_data[uid_]["fed"]
        chempot = dict(zip(fed_.chempot_diagram.elements, (x, y)))
        pts = fed_.get_transitions(chempot)
        y_shift = pts[0][1] - data_["y"][0]
        uids_to_shift = {uid_, f"{uid_}:slope"}
        _shift_slope_labels(formation_en_fig, uids_to_shift, y_shift)


@app.callback(
    Output("formation_en", "figure"),
    Input("chempot", "clickData"),
    Input("chemsys_dropdown", "value"),
    State("formation_en", "figure"),
)
def update_fed(d_click_data: dict, chemsys: str, formation_en_fig: dict) -> go.Figure:
    """Callback to update the formation energy diagram."""
    ic(ctx.triggered_id, chemsys)
    if formation_en_fig is None or ctx.triggered_id == "chemsys_dropdown":
        return _plot_fed_data(FEDS_DATA, chemsys=chemsys)
    if len(chemsys.split("-")) == 2:
        _update_fed_from_2d_chempot(
            FEDS_DATA,
            chemsys=chemsys,
            d_click_data=d_click_data,
            formation_en_fig=formation_en_fig,
        )
    return formation_en_fig


app.run(debug=True, use_reloader=False, port=8051)
# %%
