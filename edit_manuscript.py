"""Apply the value-updates that follow from the new_data (3,276-ASV) / mature-scope
re-analysis to the authoritative manuscript revision, preserving run formatting.

ONLY updates numbers whose CLAIM STRUCTURE is unchanged. The two conclusion-level
changes (p118 methanogen cross-correlation, p123 H2-correlate finding) are NOT
touched here -- they are documented in MANUSCRIPT_CHANGES.md for author review.
"""
import docx

DOCX = "manuscript_1/Co_diffusion manuscript 04.13.26.docx"


def replace_in_paragraph(para, old, new):
    """Replace `old` with `new` in a paragraph even when it is split across runs;
    the replacement inherits the first overlapping run's formatting."""
    runs = para.runs
    text = "".join(r.text for r in runs)
    idx = text.find(old)
    if idx < 0:
        return False
    end = idx + len(old)
    pos = 0
    s_run = s_off = e_run = e_off = None
    for ri, r in enumerate(runs):
        rlen = len(r.text)
        if s_run is None and pos + rlen > idx:
            s_run, s_off = ri, idx - pos
        if pos + rlen >= end:
            e_run, e_off = ri, end - pos
            break
        pos += rlen
    if s_run == e_run:
        runs[s_run].text = runs[s_run].text[:s_off] + new + runs[s_run].text[e_off:]
    else:
        runs[s_run].text = runs[s_run].text[:s_off] + new
        for ri in range(s_run + 1, e_run):
            runs[ri].text = ""
        runs[e_run].text = runs[e_run].text[e_off:]
    return True


# (paragraph_index, old, new, note)
EDITS = [
    (110, "3,216 ASVs", "3,276 ASVs", "total ASV count (new_data)"),
    (21,  "11.3% to 80.2%", "13.1% to 77.8%", "Methanobacteriaceae enrichment (highlight)"),
    (11,  "80.2% summed", "77.8% summed", "abstract: Methanobacteriaceae summed rel. abundance"),
    (111, "2.4 ± 0.3", "2.5 ± 0.2", "Shannon post-inoculum mean"),
    (111, "a low of 1.9 in", "a low of 2.1 in", "Shannon final-sample low"),
    (111, "11.3% in the inoculum to 80.2% on day 300",
          "13.1% in the inoculum to 77.8% on day 300", "Methanobacteriaceae enrichment"),
    (111, "(52% by day 300)", "(50% by day 300)", "Methanobacterium.1 day-300 abundance"),
    (129, "11.3% in the inoculum to 80.2% in the final sample",
          "13.1% in the inoculum to 77.8% in the final sample", "conclusion enrichment"),
]


def main():
    d = docx.Document(DOCX)
    P = d.paragraphs
    ok, fail = 0, 0
    for pi, old, new, note in EDITS:
        if replace_in_paragraph(P[pi], old, new):
            print(f"  [OK]   p{pi}: '{old}' -> '{new}'  ({note})")
            ok += 1
        else:
            print(f"  [FAIL] p{pi}: could not find '{old}'  ({note})")
            fail += 1
    d.save(DOCX)
    print(f"\napplied {ok} edits, {fail} failed; saved {DOCX}")


if __name__ == "__main__":
    main()
