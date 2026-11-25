# -*- coding: utf-8 -*-
"""Gera uma imagem polida do pipeline de arquitetura do Smiles Validator."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "pipeline.png")

COL_INPUT = "#37474F"   # entrada/saída
COL_OURS = "#1565C0"    # componentes do projeto
COL_BACKEND = "#6A1B9A" # backend químico
TXT = "white"

fig, ax = plt.subplots(figsize=(8.0, 11.2))
ax.set_xlim(0, 10)
ax.set_ylim(-3.4, 24.6)
ax.axis("off")

W, X = 7.0, 1.0

def box(y, title, subtitle, color, tier=None, h=2.0):
    fb = FancyBboxPatch((X, y), W, h,
                        boxstyle="round,pad=0.08,rounding_size=0.25",
                        linewidth=0, facecolor=color)
    ax.add_patch(fb)
    cx = X + W / 2
    if subtitle:
        ax.text(cx, y + h * 0.60, title, ha="center", va="center",
                color=TXT, fontsize=12, fontweight="bold")
        ax.text(cx, y + h * 0.24, subtitle, ha="center", va="center",
                color="#E8EAF6", fontsize=8.3, style="italic")
    else:
        ax.text(cx, y + h * 0.5, title, ha="center", va="center",
                color=TXT, fontsize=12, fontweight="bold")
    if tier:
        ax.text(X + W + 0.15, y + h / 2, tier, ha="left", va="center",
                fontsize=8.0, color="#37474F", fontfamily="monospace")

def arrow(y_from, y_to, label=None):
    cx = X + W / 2
    a = FancyArrowPatch((cx, y_from), (cx, y_to),
                        arrowstyle="-|>", mutation_scale=20,
                        linewidth=2.0, color="#455A64")
    ax.add_patch(a)
    if label:
        ax.text(cx + 0.2, (y_from + y_to) / 2, label, ha="left", va="center",
                fontsize=8.0, color="#37474F", style="italic")

ax.text(5.0, 23.9, "Pipeline de Validação — Smiles Validator",
        ha="center", va="center", fontsize=14.5, fontweight="bold", color="#1A237E")

box(21.9, "Cadeia SMILES", "entrada — ex.: c1ccccc1", COL_INPUT, h=1.5)
arrow(21.9, 20.1)

box(18.0, "Parser LALR(1)", "léxico com estados (SmilesLex/BracketLex) + gramática SLY",
    COL_OURS, tier="tiers: lex · grammar")
arrow(18.0, 16.2, "árvore reduzida")

box(14.1, "Semântica de anéis", "ParserManager — OpenSMILES §3.4 / §3.6",
    COL_OURS, tier="tier: ring_semantics")
arrow(14.1, 12.3)

box(10.2, "Montagem do grafo", "GraphBuilder → MolecularGraph", COL_OURS)
arrow(10.2, 8.4, "MolecularGraph")

box(6.3, "Adaptador neutro", "to_descriptor → (atoms, bonds)", COL_OURS)
arrow(6.3, 4.5, "descritor")

box(2.4, "Motor químico (backend)", "valência + aromaticidade — pysmiles / rdkit",
    COL_BACKEND, tier="tier: chemistry")
arrow(2.4, 0.6)

box(-1.0, "ValidationResult", "valid · tier · message · position   (+ tupla legada)",
    COL_INPUT, h=1.6)

leg = [(COL_OURS, "Projeto (sintaxe + grafo)"),
       (COL_BACKEND, "Backend (química)"),
       (COL_INPUT, "Entrada / Saída")]
lx = 0.6
for color, label in leg:
    ax.add_patch(FancyBboxPatch((lx, -3.1), 0.5, 0.5,
                 boxstyle="round,pad=0.02,rounding_size=0.1",
                 linewidth=0, facecolor=color))
    ax.text(lx + 0.62, -2.85, label, ha="left", va="center",
            fontsize=8.2, color="#263238")
    lx += 0.62 + len(label) * 0.135 + 0.3

plt.tight_layout()
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print("OK:", OUT)
