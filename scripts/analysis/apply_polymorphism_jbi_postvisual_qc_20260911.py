#!/usr/bin/env python3
"""Apply the final layout-only fixes found by 168-mm JBI visual inspection, fail closed.

Scientific data, statistics, axis scales and the >=9.5 pt source typography floor are
unchanged. This patch only reflows Figure 1a and removes pointwise p-value text from
Figure 3c, where the exact values remain in Results and frozen figure data.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "scripts/analysis/render_polymorphism_paper_v0_1_figures_20260910.py"

OLD_FLOW = '''    boxes = [
        (0.03, 0.62, 0.26, 0.22, "100 photos\\nper species"),
        (0.37, 0.62, 0.26, 0.22, "measurement\\n+ status gate"),
        (0.71, 0.62, 0.26, 0.22, "4 colour\\nstates"),
        (0.37, 0.16, 0.26, 0.22, "species diversity\\nD = 1 − Σp²"),
    ]
    for x, y, w, h, label in boxes:
        rect = plt.Rectangle((x, y), w, h, transform=ax.transAxes, fill=False, lw=1.1, ec=DARK)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center", fontsize=9.5)
    arrows = [((0.29, 0.73), (0.37, 0.73)), ((0.63, 0.73), (0.71, 0.73)), ((0.84, 0.62), (0.54, 0.38))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", lw=1.0, color=DARK))
    ax.text(0.5, 0.01, "Amount of polymorphism is the focal species-level estimand", transform=ax.transAxes, ha="center", fontsize=9.5)
'''

NEW_FLOW = '''    boxes = [
        (0.16, 0.78, 0.68, 0.14, "100 photos per species"),
        (0.16, 0.55, 0.68, 0.14, "measurement + status gate"),
        (0.16, 0.32, 0.68, 0.14, "4 admitted colour states"),
        (0.16, 0.07, 0.68, 0.16, "focal species diversity\\nD = 1 − Σp²"),
    ]
    for x, y, w, h, label in boxes:
        rect = plt.Rectangle((x, y), w, h, transform=ax.transAxes, fill=False, lw=1.1, ec=DARK)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, transform=ax.transAxes, ha="center", va="center", fontsize=9.5)
    arrows = [((0.50, 0.78), (0.50, 0.69)), ((0.50, 0.55), (0.50, 0.46)), ((0.50, 0.32), (0.50, 0.23))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, xycoords=ax.transAxes, arrowprops=dict(arrowstyle="->", lw=1.0, color=DARK))
'''

OLD_PTEXT = '''        for x, y, (tr, resp, _) in zip(xs, ys, groups):
            p = float(robust[(robust.tranche==tr) & (robust.response==resp) & (robust.D_variant==v)].p.iloc[0])
            ax.text(x, y + 0.006, f"{p:.3f}", ha="center", va="bottom", fontsize=9.5, color=c)
'''

NEW_PTEXT = '''        # Exact geometry-preserving P values remain in Results and the frozen figure-data table.
'''

OLD_CAPTION = "P values use the original geometry-preserving within-species spatial null arrays."
NEW_CAPTION = "Exact P values are reported in Results and retained in the frozen figure-data table; all use the original geometry-preserving within-species spatial null arrays."


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"post-JBI visual patch expected exactly one {label} target, found {count}")
    return text.replace(old, new, 1)


def main() -> int:
    text = PATH.read_text(encoding="utf-8")
    original = text
    text = replace_once(text, OLD_FLOW, NEW_FLOW, "Figure 1a flow")
    text = replace_once(text, OLD_PTEXT, NEW_PTEXT, "Figure 3c p-text loop")
    text = replace_once(text, OLD_CAPTION, NEW_CAPTION, "Figure 3 caption")

    # Keep the JBI-compliant typography/panel conventions installed by the previous patch.
    for token in [
        '"font.size": 9.5,',
        '"xtick.labelsize": 9.5,',
        '"ytick.labelsize": 9.5,',
        '"legend.fontsize": 9.5,',
        'f"({letter.lower()})"',
    ]:
        if token not in text:
            raise RuntimeError(f"JBI compliance token missing after layout patch: {token}")
    if 'fontsize=6.' in text or 'fontsize=7.' in text or 'fontsize=8.' in text:
        raise RuntimeError("sub-floor explicit font unexpectedly reintroduced")
    if text == original:
        raise RuntimeError("post-JBI visual patch made no change")
    PATH.write_text(text, encoding="utf-8", newline="\n")
    print("patched Figure 1a vertical flow and removed Figure 3c pointwise p-value labels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
