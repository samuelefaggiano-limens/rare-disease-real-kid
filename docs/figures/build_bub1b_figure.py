"""Generates docs/figures/bub1b_variants_explained.html (Step 3 figure, EN/ES).

The figure is static HTML + SVG. The language switch (top right) is pure CSS:
two hidden <input type=radio> and one rule that hides the texts with class
.es or .en depending on which one is checked. There is no JavaScript.
English is the default.

Usage:
    python docs/figures/build_bub1b_figure.py                # HTML with the EN/ES switch
    python docs/figures/build_bub1b_figure.py --static en    # English only (for PNG export)
    python docs/figures/build_bub1b_figure.py --static es    # Spanish only

The data (exon, letter triplet, protein position, reads) come from
track1_step3b_regulatory_annotated.vcf.gz; the neighbouring letters from the
GRCh38 reference genome (Ensembl). See docs/02-findings.md.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------- strings
T = {
    "page_title": ("Las dos variantes en BUB1B, explicadas", "The two BUB1B variants, explained"),
    "eyebrow": ("Track 1 · Paso 3 · gen BUB1B", "Track 1 · Step 3 · BUB1B gene"),
    "h1_a": ("Dos cambios de una sola letra en el gen ", "Two single-letter changes in the "),
    "h1_b": (", y por qué importa que haya dos", " gene, and why it matters that there are two"),
    "sub": ("El ADN es un texto escrito con cuatro letras (A, C, G, T). Un gen es un capítulo de ese texto que la célula lee de a tres letras para armar una proteína, ladrillo a ladrillo. En el genoma del paciente encontramos dos letras cambiadas en el mismo capítulo.",
            "DNA is a text written with four letters (A, C, G, T). A gene is a chapter of that text which the cell reads three letters at a time to build a protein, brick by brick. In the patient's genome we found two changed letters in the same chapter."),
    "desc": ("Mapa del gen BUB1B en el cromosoma 15 con dos variantes heterocigotas, una en el exón 17 que crea una señal de STOP prematura y otra en el exón 23 que reemplaza un aminoácido, seguido de una explicación de la heterocigosis y de los dos escenarios posibles de fase.",
             "Map of the BUB1B gene on chromosome 15 with two heterozygous variants, one in exon 17 creating a premature STOP signal and one in exon 23 replacing an amino acid, followed by an explanation of heterozygosity and the two possible phase scenarios."),
    "lang_es": ("Español", "Español"), "lang_en": ("English", "English"),

    "secA": ("A · DÓNDE ESTÁ EL GEN", "A · WHERE THE GENE IS"),
    "chrom": ("Cromosoma 15", "Chromosome 15"),
    "chrom_sub": ("1 de los 23 pares", "1 of the 23 pairs"),
    "chrom_len": ("102 MILLONES DE LETRAS", "102 MILLION LETTERS"),
    "zoom": ("ZOOM × 1700", "ZOOM × 1700"),
    "gene": ("Gen BUB1B", "BUB1B gene"),
    "gene_sub": ("60 mil letras · 23 exones", "60,000 letters · 23 exons"),
    "read_dir": ("la célula lee de izquierda a derecha, de a 3 letras", "the cell reads left to right, 3 letters at a time"),
    "change1": ("CAMBIO 1", "CHANGE 1"), "change2": ("CAMBIO 2", "CHANGE 2"),
    "aside": ("exones dibujados a igual tamaño; el gen real no es así de prolijo", "exons drawn at equal size; the real gene is not this tidy"),

    "secB": ("B · ZOOM A LAS LETRAS Y SUS VECINAS", "B · ZOOM INTO THE LETTERS AND THEIR NEIGHBOURS"),
    "p1_title": ("Exón 17 · una letra crea una señal de STOP", "Exon 17 · one letter creates a STOP signal"),
    "p1_trio": ("trío n.º 751 de 1064", "triplet no. 751 of 1064"),
    "healthy": ("SANO", "HEALTHY"), "patient": ("PACIENTE", "PATIENT"),
    "p1_h": ("= ladrillo “leucina”", "= brick “leucine”"),
    "p1_p": ("= señal de STOP: la proteína se corta en el ladrillo 751", "= STOP signal: the protein is cut short at brick 751"),
    "p1_foot": ("ClinVar: patogénico para MVA tipo 1 · en la población general: 3 de cada 100.000", "ClinVar: pathogenic for MVA type 1 · in the general population: 3 in 100,000"),
    "p2_title": ("Exón 23 (el último) · se reemplaza un ladrillo", "Exon 23 (the last) · one brick is swapped"),
    "p2_trio": ("trío n.º 1016 de 1064", "triplet no. 1016 of 1064"),
    "p2_h": ("= ladrillo “asparagina”", "= brick “asparagine”"),
    "p2_p": ("= ladrillo “lisina”: la proteína se completa con otra pieza", "= brick “lysine”: the protein is completed with a different piece"),
    "p2_foot": ("Nunca vista en la población · CADD, REVEL y AlphaGenome: alto impacto predicho", "Never seen in the population · CADD, REVEL and AlphaGenome: high predicted impact"),

    "secC": ("C · HETEROCIGOSIS: QUÉ VEMOS Y QUÉ NO PODEMOS VER", "C · HETEROZYGOSITY: WHAT WE SEE AND WHAT WE CANNOT"),
    "c_title": ("Todos tenemos dos copias de cada gen: una de cada progenitor.", "We all carry two copies of every gene: one from each parent."),
    "c_l1a": ("El secuenciador lee cada zona del ADN muchas veces. Cada cambio aparece en la ", "The sequencer reads every stretch of DNA many times. Each change shows up in "),
    "c_l1b": ("mitad", "half"),
    "c_l1c": (" de las lecturas: 25 de 46 para el cambio 1, 13 de 28 para el cambio 2.", " of the reads: 25 of 46 for change 1, 13 of 28 for change 2."),
    "c_l2a": ("Eso significa que cada cambio está en ", "That means each change sits on "),
    "c_l2b": ("una sola", "only one"),
    "c_l2c": (" de las dos copias: es ", " of the two copies: it is "),
    "c_l2d": ("heterocigoto", "heterozygous"),
    "c_l2e": (". Si estuviera en las dos, aparecería en casi todas las lecturas.", ". If it were on both, it would appear in almost every read."),
    "scA": ("ESCENARIO A · CADA CAMBIO EN UNA COPIA DISTINTA", "SCENARIO A · EACH CHANGE ON A DIFFERENT COPY"),
    "scB": ("ESCENARIO B · LOS DOS CAMBIOS EN LA MISMA COPIA", "SCENARIO B · BOTH CHANGES ON THE SAME COPY"),
    "copy1": ("Copia de un progenitor", "Copy from one parent"),
    "copy2": ("Copia del otro progenitor", "Copy from the other parent"),
    "nochange": ("SIN CAMBIOS", "NO CHANGES"),
    "scA_t": ("Ninguna copia queda intacta.", "No intact copy is left."),
    "scA_1": ("Es la “heterocigosis compuesta”: el patrón descrito en MVA tipo 1,", "This is “compound heterozygosity”: the pattern described in MVA type 1,"),
    "scA_2": ("una enfermedad que necesita que fallen las dos copias.", "a disease that requires both copies to fail."),
    "scB_t": ("Una copia sana sigue funcionando.", "One healthy copy keeps working."),
    "scB_1": ("Con el ADN de una sola persona no se puede distinguir A de B:", "With one person's DNA alone, A cannot be told apart from B:"),
    "scB_2": ("hace falta el ADN de los padres o una lectura larga que una los dos puntos.", "it takes the parents' DNA or a long read spanning both points."),

    "legend": ("LEYENDA", "LEGEND"),
    "lg_exon": ("Exón: parte del gen que se traduce a proteína", "Exon: part of the gene translated into protein"),
    "lg_intron": ("Intrón: texto intermedio que se recorta", "Intron: in-between text that is cut out"),
    "lg_change": ("Cambio encontrado en el paciente", "Change found in the patient"),
    "lg_letter": ("Letra cambiada (letras vecinas en gris)", "Changed letter (neighbouring letters in grey)"),
    "colophon": ("Fuente: VCF del paciente anotado con snpEff (GRCh38.99, transcrito canónico) + MyVariant.info (gnomAD, CADD, REVEL, ClinVar) + AlphaGenome Atlas. Letras vecinas tomadas del genoma de referencia GRCh38 (Ensembl). Anotación de investigación, no un diagnóstico clínico. Figura del Track 1, Paso 3.",
                 "Source: patient VCF annotated with snpEff (GRCh38.99, canonical transcript) + MyVariant.info (gnomAD, CADD, REVEL, ClinVar) + AlphaGenome Atlas. Neighbouring letters from the GRCh38 reference genome (Ensembl). Research annotation, not a clinical diagnosis. Track 1, Step 3 figure."),
}

MODE = "toggle"  # or "es" / "en"


def t(key: str, attrs: str = "", inner: bool = False) -> str:
    """One SVG <text> per language (or a single one in static mode)."""
    es, en = T[key]
    if MODE == "es":
        return f'<text {attrs}>{html.escape(es)}</text>'
    if MODE == "en":
        return f'<text {attrs}>{html.escape(en)}</text>'
    return (f'<text class="es" {attrs}>{html.escape(es)}</text>'
            f'<text class="en" {attrs}>{html.escape(en)}</text>')


def s(key: str) -> str:
    """Plain HTML text per language, wrapped in <span>."""
    es, en = T[key]
    if MODE == "es":
        return html.escape(es)
    if MODE == "en":
        return html.escape(en)
    return f'<span class="es">{html.escape(es)}</span><span class="en">{html.escape(en)}</span>'


MONO = "font-family=\"'Geist Mono', monospace\""
SANS = "font-family=\"'Geist', sans-serif\""
SERIF = "font-family=\"'Instrument Serif', serif\""


def letters(xs, seq, focal_idx, box_x, y=336):
    """Letter row: neighbours grey, triplet black, changed letter white on a box."""
    out = []
    for i, (x, ch) in enumerate(zip(xs, seq)):
        if i == focal_idx:
            continue
        col, w = ("#2d3142", ' font-weight="600"') if ch.isupper() and i in focal_idx_set(focal_idx) else ("#7a8399", "")
        out.append(f'<text x="{x}" y="{y}" fill="{col}"{w}>{ch}</text>')
    return "".join(out)


def focal_idx_set(i):
    return {i - 1, i, i + 1}


def zoom_panel(x0, n, seq_chars, focal_i, codon, key_prefix):
    """Panel B: 440 wide. seq_chars are 15 tokens ('…' included), 24px pitch.

    focal_i is the index of the changed letter; codon, the 3 indices of the
    triplet (the changed letter can be its 1st, 2nd or 3rd base).
    """
    xs = [x0 + 36 + 24 * i for i in range(len(seq_chars))]
    row = []
    for i, (x, ch) in enumerate(zip(xs, seq_chars)):
        if i == focal_i:
            continue
        if i in codon:
            row.append(f'<text x="{x}" y="336" fill="#2d3142" font-weight="600">{ch}</text>')
        else:
            row.append(f'<text x="{x}" y="336" fill="#7a8399">{ch}</text>')
    fx = xs[focal_i]
    bracket_l, bracket_r = xs[min(codon)] - 8, xs[max(codon)] + 8
    cx = (bracket_l + bracket_r) // 2
    return f"""
      <rect x="{x0}" y="272" width="440" height="176" rx="6" fill="#ffffff" stroke="rgba(45,49,66,0.14)" stroke-width="1"/>
      <circle cx="{x0 + 22}" cy="294" r="6" fill="#eb6c36"/>
      <text x="{x0 + 22}" y="297" fill="#ffffff" font-size="8" font-weight="600" {MONO} text-anchor="middle">{n}</text>
      {t(key_prefix + "_title", f'x="{x0 + 36}" y="298" fill="#2d3142" font-size="12" font-weight="600" {SANS}')}
      <g {MONO} font-size="12" text-anchor="middle">{"".join(row)}</g>
      <rect x="{fx - 10}" y="322" width="20" height="20" rx="3" fill="#2d3142"/>
      <text x="{fx}" y="336" fill="#ffffff" font-size="12" font-weight="600" {MONO} text-anchor="middle">{seq_chars[focal_i]}</text>
      <path d="M {bracket_l} 346 L {bracket_l} 352 L {bracket_r} 352 L {bracket_r} 346" fill="none" stroke="#2d3142" stroke-width="1"/>
      {t(key_prefix + "_trio", f'x="{cx}" y="366" fill="#4f5d75" font-size="9" {MONO} text-anchor="middle"')}
      {t("healthy", f'x="{x0 + 24}" y="396" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.14em"')}
      {t("patient", f'x="{x0 + 24}" y="416" fill="#eb6c36" font-size="8" {MONO} letter-spacing="0.14em"')}
    """


def codon_rows(x0, healthy, mutant, key_prefix):
    return f"""
      <text x="{x0 + 80}" y="396" fill="#2d3142" font-size="11" {MONO}>{healthy}</text>
      {t(key_prefix + "_h", f'x="{x0 + 128}" y="396" fill="#2d3142" font-size="11" {SANS}')}
      <text x="{x0 + 80}" y="416" fill="#2d3142" font-size="11" font-weight="600" {MONO}>{mutant}</text>
      {t(key_prefix + "_p", f'x="{x0 + 128}" y="416" fill="#2d3142" font-size="11" {SANS}')}
      {t(key_prefix + "_foot", f'x="{x0 + 24}" y="438" fill="#4f5d75" font-size="9" {SANS}')}
    """


def scenario(x0, key, copy2_markers, copy2_empty, t_key, l1, l2, m1x, m2x):
    marks = ""
    for cx, n in copy2_markers:
        marks += (f'<circle cx="{cx}" cy="665" r="5" fill="#eb6c36"/>'
                  f'<text x="{cx}" y="668" fill="#ffffff" font-size="7" font-weight="600" {MONO} text-anchor="middle">{n}</text>')
    top = ""
    for cx, n in ((m1x, 1),) if m2x is None else ((m1x, 1), (m2x, 2)):
        top += (f'<circle cx="{cx}" cy="633" r="5" fill="#eb6c36"/>'
                f'<text x="{cx}" y="636" fill="#ffffff" font-size="7" font-weight="600" {MONO} text-anchor="middle">{n}</text>')
    empty = ""
    if copy2_empty:
        empty = t("nochange", f'x="{x0 + 300}" y="669" fill="#7a8399" font-size="8" {MONO} text-anchor="middle" letter-spacing="0.14em"')
    fill2 = "#ffffff" if copy2_empty else "rgba(45,49,66,0.05)"
    return f"""
      <rect x="{x0}" y="576" width="440" height="160" rx="6" fill="#ffffff" stroke="rgba(45,49,66,0.14)" stroke-width="1"/>
      {t(key, f'x="{x0 + 20}" y="600" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.14em"')}
      {t("copy1", f'x="{x0 + 20}" y="636" fill="#4f5d75" font-size="9" {SANS}')}
      <rect x="{x0 + 180}" y="626" width="240" height="14" rx="3" fill="rgba(45,49,66,0.05)" stroke="#4f5d75" stroke-width="1"/>
      {top}
      {t("copy2", f'x="{x0 + 20}" y="668" fill="#4f5d75" font-size="9" {SANS}')}
      <rect x="{x0 + 180}" y="658" width="240" height="14" rx="3" fill="{fill2}" stroke="#4f5d75" stroke-width="1"/>
      {marks}{empty}
      {t(t_key, f'x="{x0 + 20}" y="700" fill="#2d3142" font-size="11" font-weight="600" {SANS}')}
      {t(l1, f'x="{x0 + 20}" y="714" fill="#4f5d75" font-size="9.5" {SANS}')}
      {t(l2, f'x="{x0 + 20}" y="727" fill="#4f5d75" font-size="9.5" {SANS}')}
    """


def exon_blocks():
    out = []
    for i in range(23):
        x = 188 + 32 * i
        if i in (16, 22):
            out.append(f'<rect x="{x}" y="164" width="16" height="16" rx="2" fill="rgba(235,108,54,0.08)" stroke="#eb6c36" stroke-width="1.2"/>')
        else:
            out.append(f'<rect x="{x}" y="164" width="16" height="16" rx="2" fill="#ffffff" stroke="#2d3142" stroke-width="1"/>')
    return "".join(out)


def rich_line(y, parts):
    """Section C line with bold parts: parts = [(key, bold), ...] in a single <text> per language."""
    def build(lang_i):
        spans = []
        for key, bold in parts:
            txt = html.escape(T[key][lang_i])
            spans.append(f'<tspan font-weight="600" fill="#2d3142">{txt}</tspan>' if bold else txt)
        return "".join(spans)
    attrs = f'x="40" y="{y}" fill="#4f5d75" font-size="10" {SANS}'
    if MODE == "es":
        return f"<text {attrs}>{build(0)}</text>"
    if MODE == "en":
        return f"<text {attrs}>{build(1)}</text>"
    return f'<text class="es" {attrs}>{build(0)}</text><text class="en" {attrs}>{build(1)}</text>'


def build() -> str:
    lang_attr = "es" if MODE == "es" else "en"
    switch = ""
    switch_css = ""
    if MODE == "toggle":
        switch = f"""
    <input type="radio" name="lang" id="lang-en" checked>
    <input type="radio" name="lang" id="lang-es">"""
        switch_css = """
    input[name="lang"] { position: absolute; opacity: 0; pointer-events: none; }
    #lang-en:checked ~ .frame .es, #lang-es:checked ~ .frame .en { display: none; }
    .lang-switch { position: absolute; top: 0; right: 0; display: inline-flex; border: 1px solid rgba(45,49,66,0.25);
                   border-radius: 6px; overflow: hidden; font-family: var(--font-mono); font-size: 0.66rem;
                   letter-spacing: 0.08em; text-transform: uppercase; }
    .lang-switch label { padding: 0.35rem 0.7rem; cursor: pointer; color: var(--color-muted); background: #ffffff; }
    .lang-switch label + label { border-left: 1px solid rgba(45,49,66,0.25); }
    #lang-es:checked ~ .frame .lang-switch label[for="lang-es"],
    #lang-en:checked ~ .frame .lang-switch label[for="lang-en"] { background: var(--color-ink); color: #ffffff; }
    .frame { position: relative; }"""
        switch_html = """
    <div class="lang-switch" aria-label="Language / Idioma">
      <label for="lang-en">English</label><label for="lang-es">Español</label>
    </div>"""
    else:
        switch_html = ""

    # Neighbouring letters from the GRCh38 reference genome (Ensembl), 6 on each side.
    # Change 1: ...CAGAGT[T]AAGTGC...  triplet TTA (indices 6-8), changed letter = 2nd of the triplet
    # Change 2: ...TCTGAA[T]GCCAAT...  triplet AAT (indices 5-7), changed letter = 3rd of the triplet
    p1_seq = ["…", "C", "A", "G", "A", "G", "T", "T", "A", "A", "G", "T", "G", "C", "…"]
    p2_seq = ["…", "T", "C", "T", "G", "A", "A", "T", "G", "C", "C", "A", "A", "T", "…"]

    svg = f"""
    <svg viewBox="0 0 1000 820" xmlns="http://www.w3.org/2000/svg" role="img" aria-labelledby="bub1b-title bub1b-desc">
      <title id="bub1b-title">{html.escape(T["page_title"][0 if lang_attr == "es" else 1])}</title>
      <desc id="bub1b-desc">{html.escape(T["desc"][0 if lang_attr == "es" else 1])}</desc>
      <defs>
        <marker id="arrow" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#4f5d75"/></marker>
        <marker id="arrow-accent" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#eb6c36"/></marker>
        <marker id="arrow-link" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#2e5aa8"/></marker>
      </defs>
      <rect width="100%" height="100%" fill="#f5f5f5"/>

      <!-- A. Where the gene is -->
      {t("secA", f'x="40" y="28" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.18em"')}
      {t("chrom", f'x="40" y="66" fill="#2d3142" font-size="12" font-weight="600" {SANS}')}
      {t("chrom_sub", f'x="40" y="80" fill="#4f5d75" font-size="9" {MONO}')}
      <rect x="180" y="56" width="760" height="20" rx="10" fill="rgba(45,49,66,0.05)" stroke="#4f5d75" stroke-width="1"/>
      <rect x="312" y="58" width="12" height="16" rx="6" fill="rgba(45,49,66,0.18)"/>
      <rect x="472" y="56" width="8" height="20" rx="2" fill="#2d3142"/>
      <text x="476" y="46" fill="#2d3142" font-size="9" font-weight="600" {SANS} text-anchor="middle">BUB1B</text>
      {t("chrom_len", f'x="940" y="92" fill="#7a8399" font-size="8" {MONO} text-anchor="end" letter-spacing="0.06em"')}
      <line x1="476" y1="76" x2="476" y2="112" stroke="rgba(45,49,66,0.40)" stroke-width="1" stroke-dasharray="4,3"/>
      <line x1="180" y1="112" x2="940" y2="112" stroke="rgba(45,49,66,0.40)" stroke-width="1" stroke-dasharray="4,3"/>
      <line x1="180" y1="112" x2="180" y2="148" stroke="rgba(45,49,66,0.40)" stroke-width="1" stroke-dasharray="4,3"/>
      <line x1="940" y1="112" x2="940" y2="148" stroke="rgba(45,49,66,0.40)" stroke-width="1" stroke-dasharray="4,3"/>
      <rect x="520" y="92" width="120" height="14" rx="2" fill="#f5f5f5"/>
      {t("zoom", f'x="580" y="102" fill="#4f5d75" font-size="8" {MONO} text-anchor="middle" letter-spacing="0.14em"')}
      {t("gene", f'x="40" y="176" fill="#2d3142" font-size="12" font-weight="600" {SANS}')}
      {t("gene_sub", f'x="40" y="190" fill="#4f5d75" font-size="9" {MONO}')}
      <line x1="180" y1="172" x2="940" y2="172" stroke="rgba(79,93,117,0.55)" stroke-width="1"/>
      {exon_blocks()}
      <g fill="#7a8399" font-size="8" {MONO} text-anchor="middle">
        <text x="196" y="196">1</text><text x="324" y="196">5</text><text x="484" y="196">10</text>
        <text x="644" y="196">15</text><text x="708" y="196" fill="#eb6c36" font-weight="600">17</text>
        <text x="804" y="196">20</text><text x="900" y="196" fill="#eb6c36" font-weight="600">23</text>
      </g>
      <line x1="188" y1="212" x2="300" y2="212" stroke="#4f5d75" stroke-width="1" marker-end="url(#arrow)"/>
      <rect x="308" y="205" width="200" height="14" rx="2" fill="#f5f5f5"/>
      {t("read_dir", f'x="312" y="216" fill="#4f5d75" font-size="9" {SANS}')}
      <line x1="708" y1="164" x2="708" y2="128" stroke="#eb6c36" stroke-width="1"/>
      <circle cx="708" cy="128" r="6" fill="#eb6c36"/>
      <text x="708" y="131" fill="#ffffff" font-size="8" font-weight="600" {MONO} text-anchor="middle">1</text>
      <line x1="900" y1="164" x2="900" y2="128" stroke="#eb6c36" stroke-width="1"/>
      <circle cx="900" cy="128" r="6" fill="#eb6c36"/>
      <text x="900" y="131" fill="#ffffff" font-size="8" font-weight="600" {MONO} text-anchor="middle">2</text>
      <rect x="640" y="132" width="52" height="14" rx="2" fill="#f5f5f5"/>
      {t("change1", f'x="666" y="142" fill="#eb6c36" font-size="8" {MONO} text-anchor="middle" letter-spacing="0.10em"')}
      <rect x="816" y="132" width="52" height="14" rx="2" fill="#f5f5f5"/>
      {t("change2", f'x="842" y="142" fill="#eb6c36" font-size="8" {MONO} text-anchor="middle" letter-spacing="0.10em"')}
      {t("aside", f'x="960" y="216" fill="#4f5d75" font-size="12" font-style="italic" {SERIF} text-anchor="end"')}

      <!-- B. Zoom -->
      <line x1="40" y1="236" x2="960" y2="236" stroke="rgba(45,49,66,0.10)" stroke-width="0.8"/>
      {t("secB", f'x="40" y="256" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.18em"')}
      {zoom_panel(40, 1, p1_seq, 7, (6, 7, 8), "p1")}
      {codon_rows(40, "T T A", "T G A", "p1")}
      {zoom_panel(520, 2, p2_seq, 7, (5, 6, 7), "p2")}
      {codon_rows(520, "A A T", "A A G", "p2")}

      <!-- C. Heterozygosity -->
      <line x1="40" y1="472" x2="960" y2="472" stroke="rgba(45,49,66,0.10)" stroke-width="0.8"/>
      {t("secC", f'x="40" y="492" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.18em"')}
      {t("c_title", f'x="40" y="520" fill="#2d3142" font-size="12" font-weight="600" {SANS}')}
      {rich_line(538, [("c_l1a", False), ("c_l1b", True), ("c_l1c", False)])}
      {rich_line(554, [("c_l2a", False), ("c_l2b", True), ("c_l2c", False), ("c_l2d", True), ("c_l2e", False)])}
      {scenario(40, "scA", [(446, 2)], False, "scA_t", "scA_1", "scA_2", 390, None)}
      {scenario(520, "scB", [], True, "scB_t", "scB_1", "scB_2", 870, 926)}

      <!-- Legend -->
      <line x1="40" y1="764" x2="960" y2="764" stroke="rgba(45,49,66,0.10)" stroke-width="0.8"/>
      {t("legend", f'x="40" y="780" fill="#4f5d75" font-size="8" {MONO} letter-spacing="0.18em"')}
      <rect x="40" y="792" width="16" height="12" rx="2" fill="#ffffff" stroke="#2d3142" stroke-width="1"/>
      {t("lg_exon", f'x="64" y="801" fill="#4f5d75" font-size="8.5" {SANS}')}
      <line x1="300" y1="798" x2="324" y2="798" stroke="rgba(79,93,117,0.55)" stroke-width="1"/>
      {t("lg_intron", f'x="332" y="801" fill="#4f5d75" font-size="8.5" {SANS}')}
      <circle cx="560" cy="798" r="6" fill="#eb6c36"/>
      {t("lg_change", f'x="574" y="801" fill="#4f5d75" font-size="8.5" {SANS}')}
      <rect x="760" y="791" width="14" height="14" rx="3" fill="#2d3142"/>
      {t("lg_letter", f'x="782" y="801" fill="#4f5d75" font-size="8.5" {SANS}')}
    </svg>"""

    return f"""<!DOCTYPE html>
<html lang="{lang_attr}">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{html.escape(T["page_title"][0 if lang_attr == "es" else 1])}</title>
  <link href="https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&family=Geist:wght@400;500;600&family=Geist+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    :root {{
      --color-paper: #f5f5f5; --color-ink: #2d3142; --color-muted: #4f5d75; --color-accent: #eb6c36;
      --font-sans: 'Geist', system-ui, sans-serif; --font-serif: 'Instrument Serif', serif; --font-mono: 'Geist Mono', ui-monospace, monospace;
    }}
    body {{ font-family: var(--font-sans); background: var(--color-paper); color: var(--color-ink);
           min-height: 100vh; display: flex; align-items: center; justify-content: center; padding: 3rem 2rem; }}
    .frame {{ max-width: 1200px; width: 100%; }}
    .eyebrow {{ font-family: var(--font-mono); font-size: 0.66rem; font-weight: 500; letter-spacing: 0.18em;
               text-transform: uppercase; color: var(--color-muted); margin-bottom: 0.5rem; }}
    h1 {{ font-family: var(--font-serif); font-size: clamp(1.5rem, 2.4vw + 0.75rem, 2rem); font-weight: 400;
         letter-spacing: -0.02em; line-height: 1.15; margin-bottom: 0.5rem; padding-right: 12rem; }}
    .sub {{ color: var(--color-muted); font-size: 0.9rem; max-width: 60ch; margin-bottom: 1.5rem; line-height: 1.45; }}
    svg {{ width: 100%; min-width: 900px; display: block; }}
    .colophon {{ margin-top: 1.25rem; padding-top: 0.75rem; border-top: 1px solid rgba(45,49,66,0.12);
                font-family: var(--font-mono); font-size: 0.66rem; color: var(--color-muted); line-height: 1.6; }}{switch_css}
  </style>
</head>
<body>{switch}
  <div class="frame">{switch_html}
    <p class="eyebrow">{s("eyebrow")}</p>
    <h1>{s("h1_a")}<em>BUB1B</em>{s("h1_b")}</h1>
    <p class="sub">{s("sub")}</p>
{svg}
    <p class="colophon">{s("colophon")}</p>
  </div>
</body>
</html>
"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--static", choices=["es", "en"], help="a single language, no switch")
    ap.add_argument("-o", "--output", default=str(HERE / "bub1b_variants_explained.html"))
    a = ap.parse_args()
    MODE = a.static or "toggle"
    Path(a.output).write_text(build(), encoding="utf-8")
    print("written", a.output, "mode", MODE)
