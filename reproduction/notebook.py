# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.23.0",
#     "plotly>=6.0.0",
# ]
# ///

import marimo

__generated_with = "0.23.15"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import plotly.graph_objects as go

    return go, mo


@app.cell
def _(mo):
    mo.md("""
    # Reproducing the optimal sumset–difference-set exponent

    For a finite set of integers, adding every pair and subtracting every pair
    can create collections of very different sizes. The paper builds integers
    from repeated digit patterns to make addition grow faster than subtraction.
    This notebook contains the successful-run measurements used to check that
    construction and explore a more efficient digit pattern.

    ## Verdict

    **Reproduced, with a computational extension.** Fourteen successful
    Kubernetes runs reproduced all six targeted exact claim groups and
    certified an exponent of **1.999994956618507** at depth 8,388,608.
    One cancelled sizing run and four failed runs are excluded.
    """)
    return


@app.cell
def _():
    evidence = {
        "depth": [256, 512, 1024, 1536, 2048, 4096, 6144, 8192, 16384, 24576, 32768, 65536, 98304, 131072, 262144, 393216, 524288, 1048576, 1572864, 2097152, 4194304, 6291456, 8388608],
        "exponent": [1.8184545559922474, 1.9105684643274263, 1.9563728420295896, 1.9712786555544264, 1.9786210622894687, 1.9894559960801363, 1.9930094955885311, 1.9947731614324131, 1.9974000412289172, 1.998270085503816, 1.9987039283068213, 1.9993530773736898, 1.9995689908931669, 1.9996768512484413, 1.9998385123667126, 1.9998923624453602, 1.9999192800264918, 1.9999596465146094, 1.9999730992209916, 1.9999798250179677, 1.9999899129830367, 1.9999932754336096, 1.999994956618507],
        "scaled_depth": [65536, 98304, 131072, 262144, 393216, 524288, 1048576, 1572864, 2097152, 4194304, 6291456, 8388608],
        "scaled_gap": [42.39672123786295, 42.36991923812457, 42.35575316430186, 42.33301414048765, 42.32480868525454, 42.32051347068045, 42.31369629688561, 42.31126687419601, 42.31000391906127, 42.30801559705287, 42.30731356423348, 42.30695034004748],
        "predicted_scaled_limit": 42.30574306770392,
        "gadget_labels": ["Improved base-12", "Paper base-12", "Best base-32"],
        "gadget_growth": [0.9520357417415045, 0.9686229485816499, 0.9794006899892447],
        "K": [2, 4, 8, 16, 32, 64],
        "paper_depth": [88, 132, 220, 396, 748, 1452],
        "improved_depth": [60, 90, 150, 270, 510, 990],
        "blocks": [1, 2, 4, 8, 16, 32, 64, 128, 1024, 1000000],
        "block_exponent": [1.0552229125594923, 1.0559281005082646, 1.0562776705096513, 1.0564517063842631, 1.0565385378902756, 1.0565819071412317, 1.056603580154387, 1.056614413759531, 1.0566238915778086, 1.0566252440446824],
        "block_limit": 1.0566252454310125,
    }
    return (evidence,)


@app.cell
def _(evidence, go, mo):
    primary_gap = [2.0 - value for value in evidence["exponent"]]
    primary_figure = go.Figure(
        go.Scatter(
            x=evidence["depth"],
            y=primary_gap,
            mode="lines+markers",
            line={"color": "#2463A8", "width": 3},
        )
    )
    primary_figure.update_layout(
        title="Certified exponent closes in on 2",
        xaxis_title="Construction depth",
        yaxis_title="Remaining gap, 2 − certified exponent",
        template="plotly_white",
    )
    primary_figure.update_xaxes(type="log")
    primary_figure.update_yaxes(type="log")
    mo.vstack(
        [
            mo.md(
                r"""
                ## Primary result

                Lower is better: the vertical axis shows the distance remaining
                to exponent 2. The gap falls from about 0.182 at depth 256 to
                \(5.0434\\times10^{-6}\) at depth 8,388,608.
                """
            ),
            primary_figure,
        ]
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Setup and exact checks

    For an integer set \(A\), the sumset \(A+A\) contains all pairwise sums
    and the difference set \(A-A\) all pairwise differences. The measured
    exponent is

    \[
    C=\frac{\log |A+A|}{\log |A-A|}.
    \]

    The baseline passed all six targeted exact claim groups. It reproduced
    the paper digit set \(\{0,1,2,4,5,9\}\), its three-state carry transition,
    and recurrence values \(1,11,127,1475,17143\). Direct enumeration of the
    toy construction matched \(|R|=432\), modular sizes 3,600 and 3,328,
    \(|A|=864\), \(|A+A|=13,859\), and \(|A-A|=13,125\).
    """)
    return


@app.cell
def _(evidence, go, mo):
    gadget_figure = go.Figure(
        go.Bar(
            x=evidence["gadget_growth"],
            y=evidence["gadget_labels"],
            orientation="h",
            marker_color=["#18864B", "#2463A8", "#667085"],
            text=[f"{value:.4f}" for value in evidence["gadget_growth"]],
            textposition="outside",
        )
    )
    gadget_figure.update_layout(
        title="Improved base-12 gadget has the best carry behavior",
        xaxis_title="Normalized carry growth (lower is better)",
        template="plotly_white",
    )
    gadget_figure.update_xaxes(range=[0.945, 0.985])
    gadget_figure.update_yaxes(autorange="reversed")
    mo.vstack(
        [
            mo.md(
                r"""
                ## A more efficient digit gadget

                The search found digits \(\{0,1,3,4,5,8\}\), with normalized
                carry growth 0.9520357 versus the paper gadget's 0.9686229.
                An exhaustive base-32 search checked 300,540,195 normalized masks
                and carry-audited all 64 finalists; its best rate was worse.
                """
            ),
            gadget_figure,
        ]
    )
    return


@app.cell
def _(evidence, go, mo):
    efficiency_figure = go.Figure()
    efficiency_figure.add_trace(
        go.Scatter(
            x=evidence["K"],
            y=evidence["paper_depth"],
            mode="lines+markers",
            name="Paper gadget",
        )
    )
    efficiency_figure.add_trace(
        go.Scatter(
            x=evidence["K"],
            y=evidence["improved_depth"],
            mode="lines+markers",
            name="Improved gadget",
        )
    )
    efficiency_figure.update_layout(
        title="Improved carry control cuts sufficient depth by 31.8%",
        xaxis_title="Theorem parameter K",
        yaxis_title="Required construction depth",
        template="plotly_white",
    )
    efficiency_figure.update_xaxes(
        type="log",
        tickmode="array",
        tickvals=evidence["K"],
        ticktext=[str(value) for value in evidence["K"]],
    )
    mo.vstack(
        [
            efficiency_figure,
            mo.md(
                """
                The sufficient depth changes from 22(K + 2) to 15(K + 2).
                This improves efficiency without changing the limiting theorem.
                """
            ),
        ]
    )
    return


@app.cell
def _(evidence, go, mo):
    scaled_figure = go.Figure()
    scaled_figure.add_trace(
        go.Scatter(
            x=evidence["scaled_depth"],
            y=evidence["scaled_gap"],
            mode="lines+markers",
            name="Measured d(2 − C)",
        )
    )
    scaled_figure.add_hline(
        y=evidence["predicted_scaled_limit"],
        line_dash="dash",
        line_color="#D97706",
        annotation_text="Predicted limit",
    )
    scaled_figure.update_layout(
        title="Scaled gap converges to the predicted constant",
        xaxis_title="Construction depth d",
        yaxis_title="Scaled gap d(2 − C)",
        template="plotly_white",
    )
    scaled_figure.update_xaxes(type="log")
    mo.vstack(
        [
            mo.md(
                """
                ## Asymptotic diagnostic

                The scaled gap fell from 42.3967 at depth 65,536 to 42.3069503
                at depth 8,388,608, compared with the predicted limit
                42.3057431. The final residual is 0.0012073.
                """
            ),
            scaled_figure,
        ]
    )
    return


@app.cell
def _(evidence, go, mo):
    block_figure = go.Figure(
        go.Scatter(
            x=evidence["blocks"],
            y=evidence["block_exponent"],
            mode="lines+markers",
            name="Measured exponent",
        )
    )
    block_figure.add_hline(
        y=evidence["block_limit"],
        line_dash="dash",
        line_color="#D97706",
        annotation_text="Computed limit",
    )
    block_figure.update_layout(
        title="Repeated blocks approach a stable limiting exponent",
        xaxis_title="Number of repeated blocks",
        yaxis_title="Certified exponent",
        template="plotly_white",
    )
    block_figure.update_xaxes(type="log")
    mo.vstack(
        [
            mo.md(
                """
                ## Block-lift robustness

                All eight tested parameter grids were monotone. In the displayed
                case, the exponent rises from 1.0552229 at one block to
                1.0566252440 at one million blocks, matching the computed limit
                1.0566252454.
                """
            ),
            block_figure,
        ]
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ## Evidence boundary

    Every number used above came from terminal logs of 14 successful
    Kubernetes runs. A cancelled sizing run, an initial setup failure, and
    three deadline-limited larger-base searches contribute no claims. The
    finite computations certify the displayed candidates and exponents;
    they complement rather than replace the paper's general proof.
    """)
    return


if __name__ == "__main__":
    app.run()
